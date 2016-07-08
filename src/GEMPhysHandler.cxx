#include "GEMPhysHandler.h"

#include "GEMOnlineHitDecoder.h"
#include "GEMZeroHitDecoder.h"

#include "PRadDataHandler.h"
#include "PRadReconstructor.h"
#include "PRadTDCGroup.h"
#include "PRadBenchMark.h"

using namespace std;

#define PI 3.1415926

GEMPhysHandler::GEMPhysHandler()
{
  pHandler = new PRadDataHandler();
  pDSTParser = new PRadDSTParser(pHandler);

  pHandler -> ReadConfig("config.txt");

  pGEMReconstructor = new PRadGEMReconstructor(pHandler);

  totalEnergyDeposit = 1800; //MeV
  beamEnergy = 2.147;//GeV

  //compute intersection points
  px1 = 0x270F;py1=0x270F;
  px2=0x270F;py2=0x270F;
  cx1 = 0x270F;cy1=0x270F;
  cx2=0x270F;cy2=0x270F;

  px1_c = 0x270F;py1_c=0x270F;
  px2_c=0x270F;py2_c=0x270F;
  cx1_c = 0x270F;cy1_c=0x270F;
  cx2_c=0x270F;cy2_c=0x270F;
 
  //test
  hhTimeCorrelation = new TH2F("hhTimeCorrelation", "HyCal vs Scin", 1000, 0, 10000, 1000, 0, 10000);
  hTimeDiff = new TH1F("hTimeDiff", "Scin - HyCal", 1000, -5000, 5000);

  BookHistos();
}

GEMPhysHandler::~GEMPhysHandler()
{
}

void GEMPhysHandler::ProcessAllFiles()
{
  int nFile = config.nFile;

  cout<<"GEMPhysHandler::ProcessAllFiles: # of Files to analyze: "
      <<nFile
      <<endl;

  for(int i=0;i<nFile;i++)
    {
      cout<<config.fileList[i]<<endl;
    }

  for(int i=0;i<nFile;i++)
    {
      filename = config.fileList[i];
      ProcessAllEvents(-1);
      //ProcessAllEvents(40000);
    }

  // Save Histos
  SavePhysResults();
}

int GEMPhysHandler::ProcessAllEvents(int evtID )
{
  cout<<"Process File:  "<<filename<<endl;
  int entry = 0;

  PRadBenchMark timer;

  while(pDSTParser->Read() && entry < 30000)
    {
      if( pDSTParser->EventType() == PRad_DST_Event) 
	{
	  ++entry;
	  Analyze();
	}
      else if (pDSTParser ->EventType() == PRad_DST_Epics)
	{
	  pHandler -> GetEPICSData().push_back(pDSTParser->GetEPICSEvent());
	}
    }

  cout << "used time: " << timer.GetElapsedTime() << " ms" << endl;

  return entry;
}

void GEMPhysHandler::Analyzer()
{
  auto event = pDSTParser -> GetEvent();
  auto gem_cluster = pGEMReconstrctor->CoarseGEMReconstruct(event);
  auto hycal_cluster = pHandler->GetHyCalCluster(event);

  hHyCalClusterMul->Fill(hycal_cluster->size());
  for(int i=0;i<hycal_cluster->size(); i++)
    {
      hhHyCalClusterMap->Fill(hycal_cluster->at(i).x, hycal_cluster->at(i).y);
      hHyCalEnergy->Fill(hycal_cluster->at(i).E);

      if(hycal_cluster->size() == 1) 
	hHyCalEnergyEp->Fill(hycal_cluster->at(i).E);
      else if(hycal_cluster->size() == 2) 
	{
	  hHyCalEnergyMoller->Fill(hycal_cluster->at(i).E);
	}
    }

  vector<PRadGEMCluster> gem1, gem2;
  vector<PRadGEMCluster> hgem1, hgem2;
  pGEMClusterReconstruct->GEMClusteringLocal(gem1, gem2);
  pGEMClusterReconstruct->GEMClusteringHycal(hgem1, hgem2);

  online_hit.FillHistos(hNbClusterPerPlaneX, hNbClusterPerPlaneY, hClusterDistX, hClusterDistY);

  //general cut
  //if( (gem1.size()>=5 ) || (gem2.size()>=5) ) continue;
  //if( energy < 900  || energy > 1300) continue;

  // offset between GEMs
  // [1]: refer to "GEMZeroHitDecoder.cxx"
  /*
   * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   * overlapping area: 44mm
   * x side length: 550.4
   * overlapping area starting with: 550.4/2 -44 = 231.2
   *
   * origin trasfer to beam hole center:
   *    transfered distance: 550.4/2 - 44/2 = 253.2
   *    GEM1 coordinate transfer: x1 = x1 - 253.2; y1 = y1
   *    GEM2 coordinate transfer: x2 = 253.2 - x2; y2 = -y2
   * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   */

  // compute offsets between two gem chambers
  if( (gem1.size()==1) && (gem2.size()==1) )
    {
      if( (gem1[0].x > OverlapStart) && (gem2[0].x > OverlapStart) )
	{
	  {
	    //z offset correction
	    double z_gem1 = 5300; //mm
	    double z_gem2 = 5260; //mm
	    double x1prime = (gem1[0].x-O_Transfer)*z_gem2/z_gem1;
	    double y1prime = gem1[0].y*z_gem2/z_gem1;
	    hYOffsetGEM->Fill( ( y1prime ) - ( -gem2[0].y) );
	    hXOffsetGEM->Fill(( ( x1prime) - ( O_Transfer-gem2[0].x) ));

	    //after correction
	    double xoffset = -0.3618;// x1 - x2
	    double yoffset = 0.1792; // y1 - y2
	    x1prime = (gem1[0].x - xoffset -O_Transfer)*z_gem2/z_gem1;
	    y1prime = (gem1[0].y - yoffset )*z_gem2/z_gem1;
	    hYOffsetGEMAfterCorrection->Fill( ( y1prime ) - ( -gem2[0].y) );
	    hXOffsetGEMAfterCorrection->Fill(( ( x1prime) - ( O_Transfer-gem2[0].x) ));
	  }
	}
    }

  CharactorizeGEM(&online_hit);
  ProcessEp(&online_hit);
  ProcessMoller(&online_hit);
  ProcessMollerAfterCorrection(&online_hit);

  if(hgem1.size() > 0)
    {
      for(int i=0;i<hgem1.size();i++)
	{
	  hhGEMClusterMap->Fill(hgem1[i].x, hgem1[i].y);
	}
      if( hgem1.size() == 1)
	{
	  hhGEMClusterMapSingleCluster->Fill(hgem1[0].x, hgem1[0].y);
	}
      if( hgem1.size() <= 2)
	{
	  hhGEMClusterMapMaxTwoCluster->Fill(hgem1[0].x, hgem1[0].y);
	}
    }
  if(hgem2.size() > 0)
    {
      for(int i=0;i<hgem2.size();i++)
	{
	  hhGEMClusterMap->Fill(hgem2[i].x,hgem2[i].y);
	}
      if( hgem2.size() == 1)
	{
	  hhGEMClusterMapSingleCluster->Fill(hgem2[0].x, hgem2[0].y);
	}
      if( hgem2.size() <= 2)
	{
	  hhGEMClusterMapMaxTwoCluster->Fill(hgem2[0].x, hgem2[0].y);
	}
    }

}

template<class T> void GEMPhysHandler::ProcessEp(T * hit_decoder)
{
  // theta distribution
  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm
  double theta = 0;
  double q_square = 0;
  double top = 0.0;
  double bottom = 0.0; //q_square = top/bottom
  double e0 = beamEnergy; //GeV

  vector<PRadGEMCluster> gem1, gem2;
  hit_decoder->GetClusterHyCalPlusMode(gem1, gem2);
  
  int s = HyCalGEMPosMatch(gem1, gem2, hycal_cluster);
  // offset correction to beam line
  double xoffset_beam = 1.631;
  double yoffset_beam = 0.366;
  
  if( s==0 ) return;

  //ep events
  //theta distribution for single cluster
  if( (gem1.size() == 1)&&(gem2.size() == 0)) 
    { 
      // Ep energy cut
      if(gem1[0].energy < totalEnergyDeposit ) return;
      gem1[0].x -= xoffset_beam;
      gem1[0].y -= yoffset_beam;

      theta = TMath::Sqrt((gem1[0].x)*(gem1[0].x) + gem1[0].y*gem1[0].y) / z_gem1;
      theta = TMath::ATan(theta);
      double theta_d = theta*180.0/PI;
      hThetaDistributionEp->Fill(theta_d);
      hhEnergyVsAngle->Fill(theta_d, gem1[0].energy);
      hhEnergyVsAngleEp->Fill(theta_d, gem1[0].energy);

      //q suqare
      top = 4.0*beamEnergy*beamEnergy*TMath::Sin(theta/2)*TMath::Sin(theta/2);
      bottom = 1+ (2*beamEnergy/0.938)*TMath::Sin(theta/2)*TMath::Sin(theta/2);
      q_square = top/bottom;
      hQSquareEp->Fill(q_square);
      hhQSquareScattAngleEp -> Fill(theta_d, q_square);
      if( (theta_d > 0.7) && (theta_d<0.8) ) hQSquareEp1->Fill(q_square);
      if( (theta_d > 1.0) && (theta_d<1.1) ) hQSquareEp2->Fill(q_square);
      if( (theta_d > 1.5) && (theta_d<1.6) ) hQSquareEp3->Fill(q_square);
      if( (theta_d > 2.0) && (theta_d<2.1) ) hQSquareEp4->Fill(q_square);
      if(theta_d > 2.2) hQSquareEp5->Fill(q_square);

    }
  else if( (gem2.size() == 1)&&(gem1.size() == 0))
    {
      //Ep energy cut
      if(gem2[0].energy < totalEnergyDeposit) return;

      gem2[0].x -= xoffset_beam;
      gem2[0].y -= yoffset_beam;

      theta = TMath::Sqrt( (gem2[0].x)*(gem2[0].x) + gem2[0].y*gem2[0].y) / z_gem2;
      theta = TMath::ATan(theta);
      double theta_d = theta*180.0/PI;
      hThetaDistributionEp->Fill(theta_d);
      hhEnergyVsAngle->Fill(theta_d, gem2[0].energy);
      hhEnergyVsAngleEp->Fill(theta_d, gem2[0].energy);

      //q suqare
      top = 4.0*beamEnergy*beamEnergy*TMath::Sin(theta/2)*TMath::Sin(theta/2);
      bottom = 1+ (2*beamEnergy/0.938)*TMath::Sin(theta/2)*TMath::Sin(theta/2);
      q_square = top/bottom;
      hQSquareEp->Fill(q_square); 
      hhQSquareScattAngleEp -> Fill(theta_d, q_square);
      if((theta_d > 0.7) && (theta_d<0.8) ) hQSquareEp1->Fill(q_square);
      if((theta_d > 1.0)&& (theta_d<1.1) ) hQSquareEp2->Fill(q_square);
      if((theta_d > 1.5)&& (theta_d<1.6) ) hQSquareEp3->Fill(q_square);
      if((theta_d > 2.0)&& (theta_d<2.1) ) hQSquareEp4->Fill(q_square);
      if(theta_d > 2.2) hQSquareEp5->Fill(q_square);

    }
  else if ( (gem2.size() == 1)&&(gem1.size() == 1) )
    {
      /*
       * Note:
       *     coordinates returned from class ZeroHitDecoder or OnlineHitDecoder
       *     is in chamber coordinate. 
       *     Total length in X side: 550.4mm
       *
       *     Need to find out the overlap region. 51.2mm overlapping. 
       *     checked the gerber file. only the hole area is overlapped.
       *     overlapping x region: x >  550.4/2  - 51.2 (=OverlapStart)
       */
      //if( (gem1[0]>OverlapStart) && (gem2[0].x>OverlapStart)  && ((TMath::Abs(gem2[0].y+gem1[0].y))<2.0 /* suppose the shift between two chambers is smaller than 2mm  */) )
      /*
	{
	theta = TMath::Sqrt( (gem2[0].x)*(gem2[0].x) + gem2[0].y*gem2[0].y) / z_gem2;
	theta = TMath::ATan(theta);
	hThetaDistributionEp->Fill(theta*180.0/PI);

	//q suqare
	top = 4.0*1.1*1.1*TMath::Sin(theta/2)*TMath::Sin(theta/2);
	bottom = 1+ (2*1.1/0.938)*TMath::Sin(theta/2)*TMath::Sin(theta/2);
	q_square = top/bottom;
	hQSquareEp->Fill(q_square);

	}
      */
    }

}

template<class T> void GEMPhysHandler::ProcessMoller(T * hit_decoder)
{
  // theta distribution
  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm
  double theta = 0;
  double q_square = 0;
  double top = 0.0;
  double bottom = 0.0; //q_square = top/bottom
  double e0 = beamEnergy; //GeV
  
  vector<PRadGEMCluster> gem1, gem2;
  hit_decoder->GetClusterHyCalPlusMode(gem1, gem2);

  GeometryMollerRing(gem1, gem2);
 
  int s = HyCalGEMPosMatch(gem1, gem2, hycal_cluster);

  if( s == 0 ) return ;
  //Moller events
  //theta distribution for moller

  // moller 2d ring
  float ThetaSmall = 0.7/180.0 * PI;
  float ThetaLarge = 0.8/180.0 * PI;
  float ThetaSmall2 = 1.0/180.0 * PI;
  float ThetaLarge2 = 1.1/180.0 * PI;

  // origin offset from beam line
  // z coordinate: GEM2 : 5260mm
  double xoffset_beam = 1.631;
  double yoffset_beam = 0.366;


  // Moller events selection
  //     require them in different quadrants, very rough
  if( (gem1.size() == 2)&&(gem2.size() == 0)) 
    { 
      if(  ( (gem1[0].x)*(gem1[1].x) < 0 ) && ( gem1[0].y*gem1[1].y < 0 ) ) 
	{
	  if( (gem1[0].energy+gem1[1].energy) < totalEnergyDeposit) return;
	  double temp = 0;
	  double slope1 = 0;
	  double slope2 = 0;
	  //1st electron
	  theta = TMath::Sqrt((gem1[0].x)*(gem1[0].x) + gem1[0].y*gem1[0].y) / z_gem1;
	  theta = TMath::ATan(theta);
	  hhEnergyVsAngle->Fill(theta*180.0/PI, gem1[0].energy);
	  hhEnergyVsAngleMoller->Fill(theta*180.0/PI, gem1[0].energy);
	  
	  slope1 = TMath::Abs( gem1[0].y/gem1[0].x );
	  slope1 = TMath::ATan(slope1); assert( slope1 < PI/2);
	  if( (gem1[0].x<0) && (gem1[0].y>0) ) slope1 = slope1+ PI/2; 
	  else if( (gem1[0].x<0) && (gem1[0].y<0) ) slope1 = slope1+PI; 
	  else if( (gem1[0].x>0) && (gem1[0].y <0)) slope1 = slope1 + 1.5*PI;
	  //hThetaDistributionMoller->Fill(theta*180.0/PI);
	  temp = theta;
	  //2nd electron
	  theta = TMath::Sqrt((gem1[1].x)*(gem1[1].x) + gem1[1].y*gem1[1].y) / z_gem1;
	  theta = TMath::ATan(theta);
	  hhEnergyVsAngle->Fill(theta*180.0/PI, gem1[1].energy);
	  hhEnergyVsAngleMoller->Fill(theta*180.0/PI, gem1[1].energy);

	  slope2 = gem1[1].y / (gem1[1].x);
	  slope2 = TMath::ATan(slope2);
	  if( (gem1[1].x<0) && (gem1[1].y>0) ) slope1 = slope1+ PI/2; 
	  else if( (gem1[1].x<0) && (gem1[1].y<0) ) slope1 = slope1+PI; 
	  else if( (gem1[1].x>0) && (gem1[1].y <0)) slope1 = slope1 + 1.5*PI;
	  if(theta > temp) 
	    {
	      hThetaDistributionMollerLarge->Fill(theta*180.0/PI); 
	      hThetaDistributionMollerSmall->Fill(temp*180.0/PI);
	      if( ((temp*180.0/PI)>0.52) && ((temp*180.0/PI)<0.62) ) hThetaDistributionMollerSlice->Fill( theta*180.0/PI);
	    }
	  else 
	    {
	      hThetaDistributionMollerLarge->Fill(temp*180.0/PI); 
	      hThetaDistributionMollerSmall->Fill(theta*180.0/PI);
	      if( ((theta*180.0/PI)>0.54) && ((theta*180.0/PI)<0.56) ) hThetaDistributionMollerSlice->Fill( temp*180.0/PI);
	    }

	  //2d ring
	  if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge)   ) )
	    {
	      hhMoller2DRing->Fill(gem1[0].x, gem1[0].y);
	      hhMoller2DRing->Fill(gem1[1].x, gem1[1].y);
	    }
	  if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2)   ) )
	    {
	      hhMoller2DRing2->Fill(gem1[0].x, gem1[0].y);
	      hhMoller2DRing2->Fill(gem1[1].x, gem1[1].y);
	    }

	  theta+=temp;
	  hThetaDistributionMoller->Fill(theta*180.0/PI);
	  hLinearDeviationMoller->Fill((slope1-slope2)*180.0/PI - 180.0);

	  //x, y offset
	  //if(TMath::Abs( (slope1-slope2)*180.0/PI) < 1.0)
	  {
	    if(px1 != 0x270F)
	      {
	        px1 = cx1;
		py1 = cy1;
		px2 = cx2;
		py2 = cy2;
		cx1 = gem1[0].x;
		cx2 = gem1[1].y;
		cy1 = gem1[0].y;
		cy2 = gem1[1].y;

		double a1 = py1 - py2;
		double b1 = px2 - px1;
		double c1 = px1 * py2 - px2 * py1;

		double a2 = cy1 - cy2;
		double b2 = cx2 - cx1;
		double c2 = cx1 * cy2 - cx2 * cy1;
	        
		double D = a1*b2 - a2*b1;

		if(D != 0)
		  {
		    double xi = (b1*c2 - b2*c1)/D; //cout<<xi<<endl;
		    double yi = (a2*c1 - a1*c2)/D; //cout<<yi<<endl;
		    hXOffsetFromMoller->Fill(xi);
		    hYOffsetFromMoller->Fill(yi);
		    hhMollerCenter->Fill(xi, yi);
		  }

	      }
	    else
	      {
	        px1 = cx1 = gem1[0].x;
		py1 = cy1 = gem1[0].y;
		px2 = cx2 = gem1[1].x;
		py2 = cy2 = gem1[1].y;

	      }
	  }
	} 
    }
  else if( (gem2.size() == 2)&&(gem1.size() == 0))
    { 
      if(  ( (gem2[0].x)*(gem2[1].x) < 0 ) && ( gem2[0].y*gem2[1].y < 0 ) ) 
	{
	  if( (gem2[0].energy+gem2[1].energy) < totalEnergyDeposit) return;
	  double temp = 0;
	  double slope1 = 0;
	  double slope2 = 0;
	  //1st electron
	  theta = TMath::Sqrt( (gem2[0].x)*(gem2[0].x) + gem2[0].y*gem2[0].y) / z_gem2;
	  theta = TMath::ATan(theta);
	  hhEnergyVsAngle->Fill(theta*180.0/PI, gem2[0].energy);
	  hhEnergyVsAngleMoller->Fill(theta*180.0/PI, gem2[0].energy);

          slope1 = TMath::Abs( gem2[0].y/gem2[0].x );
	  slope1 = TMath::ATan(slope1); assert( slope1 < PI/2);
	  if( (gem2[0].x<0) && (gem2[0].y>0) ) slope1 = slope1+ PI/2; 
	  else if( (gem2[0].x<0) && (gem2[0].y<0) ) slope1 = slope1+PI; 
	  else if( (gem2[0].x>0) && (gem2[0].y <0)) slope1 = slope1 + 1.5*PI;
	  //hThetaDistributionMoller->Fill(theta*180.0/PI);
	  temp = theta;
	  //2nd electron
	  theta = TMath::Sqrt( (gem2[1].x)*(gem2[1].x) + gem2[1].y*gem2[1].y) / z_gem2;
	  theta = TMath::ATan(theta);
	  hhEnergyVsAngle->Fill(theta*180.0/PI, gem2[1].energy);
	  hhEnergyVsAngleMoller->Fill(theta*180.0/PI, gem2[1].energy);

	  slope2 = gem2[1].y / (gem2[1].x);
	  slope2 = TMath::ATan(slope2);
	  if( (gem2[1].x<0) && (gem2[1].y>0) ) slope1 = slope1+ PI/2; 
	  else if( (gem2[1].x<0) && (gem2[1].y<0) ) slope1 = slope1+PI; 
	  else if( (gem2[1].x>0) && (gem2[1].y <0)) slope1 = slope1 + 1.5*PI;
	  if(theta > temp) 
	    {
	      hThetaDistributionMollerLarge->Fill(theta*180.0/PI); 
	      hThetaDistributionMollerSmall->Fill(temp*180.0/PI);
	      if( ((temp*180.0/PI)>0.52) && ((temp*180.0/PI)<0.62) ) hThetaDistributionMollerSlice->Fill( theta*180.0/PI);
	    }
	  else 
	    {
	      hThetaDistributionMollerLarge->Fill(temp*180.0/PI); 
	      hThetaDistributionMollerSmall->Fill(theta*180.0/PI);
	      if( ((theta*180.0/PI)>0.54) && ((theta*180.0/PI)<0.56) ) hThetaDistributionMollerSlice->Fill( temp*180.0/PI);
	    }
	  //2d ring
	  if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge)   ) )
	    {
	      hhMoller2DRing->Fill(gem2[0].x, gem2[0].y);
	      hhMoller2DRing->Fill(gem2[1].x, gem2[1].y);
	    }
	  if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2)   ) )
	    {
	      hhMoller2DRing2->Fill(gem2[0].x, gem2[0].y);
	      hhMoller2DRing2->Fill(gem2[1].x, gem2[1].y);
	    }
	  theta+=temp;
	  hThetaDistributionMoller->Fill(theta*180.0/PI);
	  hLinearDeviationMoller->Fill((slope1-slope2)*180.0/PI - 180.0);

	  //x, y offset
	  //if(TMath::Abs( (slope1-slope2)*180.0/PI) < 1.0)
	  {
	    if(px1 != 0x270F)
	      {
	        px1 = cx1;
		py1 = cy1;
		px2 = cx2;
		py2 = cy2;
		cx1 = gem2[0].x;
		cx2 = gem2[1].x;
		cy1 = gem2[0].y;
		cy2 = gem2[1].y;

		double a1 = py1 - py2;
		double b1 = px2 - px1;
		double c1 = px1 * py2 - px2 * py1;

		double a2 = cy1 - cy2;
		double b2 = cx2 - cx1;
		double c2 = cx1 * cy2 - cx2 * cy1;
	        
		double D = a1*b2 - a2*b1;

		if(D != 0)
		  {
		    double xi = (b1*c2 - b2*c1)/D; //cout<<xi<<endl;
		    double yi = (a2*c1 - a1*c2)/D; //cout<<yi<<endl;
		    hXOffsetFromMoller->Fill(xi);
		    hYOffsetFromMoller->Fill(yi);
		    hhMollerCenter->Fill(xi, yi);
		  }

	      }
	    else
	      {
	        px1 = cx1 = gem2[0].x;
		py1 = cy1 = gem2[0].y;
		px2 = cx2 = gem2[1].x;
		py2 = cy2 = gem2[1].y;

	      }
	  }
	}
    }
  else if ( (gem2.size() == 1)&&(gem1.size() == 1) )
    {
      //if( ( TMath::Abs(gem1[0].x+gem2[0].x)>2.0)  && ((TMath::Abs(gem2[0].y+gem1[0].y))>2.0 /* suppose the shift between two chambers is smaller than 2mm  */) )
      {
	if(  ( (gem1[0].x)*(gem2[0].x) < 0 ) && ( gem1[0].y* (gem2[0].y) < 0 ) ) 
	  { 
	    if( (gem1[0].energy + gem2[0].energy) < totalEnergyDeposit) return; 
	    double temp = 0;
	    double slope1 = 0;
	    double slope2 = 0;
	    double symmetricTheta1 = 0;
	    double symmetricTheta2 = 0;
	    //1st electron
	    theta = TMath::Sqrt( (gem1[0].x)*(gem1[0].x) + gem1[0].y*gem1[0].y) / z_gem1;
	    theta = TMath::ATan(theta);
	    symmetricTheta1 = theta * 180.0/PI;
	    hhEnergyVsAngle->Fill(theta*180.0/PI, gem1[0].energy);
	    hhEnergyVsAngleMoller->Fill(theta*180.0/PI, gem1[0].energy);

	    slope1 = gem1[0].y / (gem1[0].x);
	    slope1 = TMath::ATan(slope1);
	    if( (gem1[0].x<0) && (gem1[0].y>0) ) slope1 = slope1+ PI/2; 
	    else if( (gem1[0].x<0) && (gem1[0].y<0) ) slope1 = slope1+PI; 
	    else if( (gem1[0].x>0) && (gem1[0].y <0)) slope1 = slope1 + 1.5*PI;
	    //hThetaDistributionMoller->Fill(theta*180.0/PI);
	    temp = theta;
	    //2nd electron
	    theta = TMath::Sqrt( (gem2[0].x)*(gem2[0].x) + gem2[0].y*gem2[0].y) / z_gem2;
	    theta = TMath::ATan(theta);
	    symmetricTheta2 = theta*180.0/PI;
	    hhEnergyVsAngle->Fill(theta*180.0/PI, gem2[0].energy);
	    hhEnergyVsAngleMoller->Fill(theta*180.0/PI, gem2[0].energy);

	    slope2 = gem2[0].y / (gem2[0].x);
	    slope2 = TMath::ATan(slope2);
	    if( (gem2[0].x<0) && (gem2[0].y>0) ) slope1 = slope1+ PI/2; 
	    else if( (gem2[0].x<0) && (gem2[0].y<0) ) slope1 = slope1+PI; 
	    else if( (gem2[0].x>0) && (gem2[0].y <0)) slope1 = slope1 + 1.5*PI;
	    if(theta > temp) 
	      {
		hThetaDistributionMollerLarge->Fill(theta*180.0/PI); 
		hThetaDistributionMollerSmall->Fill(temp*180.0/PI);
		if( ((temp*180.0/PI)>0.52) && ((temp*180.0/PI)<0.62) ) hThetaDistributionMollerSlice->Fill( theta*180.0/PI);
	      }
	    else 
	      {
		hThetaDistributionMollerLarge->Fill(temp*180.0/PI); 
		hThetaDistributionMollerSmall->Fill(theta*180.0/PI);
		if( ((theta*180.0/PI)>0.54) && ((theta*180.0/PI)<0.56) ) hThetaDistributionMollerSlice->Fill( temp*180.0/PI);
	      }
	    //2d ring
	    if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge)   ))
              {
	        hhMoller2DRing->Fill(gem2[0].x, gem2[0].y);
	        hhMoller2DRing->Fill(gem1[0].x, gem1[0].y);
	      }
	    if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2)   ))
              {
	        hhMoller2DRing2->Fill(gem2[0].x, gem2[0].y);
	        hhMoller2DRing2->Fill(gem1[0].x, gem1[0].y);
	      }
	    theta+=temp;
	    hThetaDistributionMoller->Fill(theta*180.0/PI);
	    hLinearDeviationMoller->Fill((slope1-slope2)*180.0/PI - 180.0);

	    //x, y offset
	    //if(TMath::Abs( (slope1-slope2)*180.0/PI) < 1.0)
	    {
	      if(px1 != 0x270F)
		{
		  px1 = cx1;
		  py1 = cy1;
		  px2 = cx2;
		  py2 = cy2;
		  cx1 = gem1[0].x;
		  cx2 = gem2[0].x;
		  cy1 = gem1[0].y;
		  cy2 = gem2[0].y;

		  double a1 = py1 - py2;
		  double b1 = px2 - px1;
		  double c1 = px1 * py2 - px2 * py1;

		  double a2 = cy1 - cy2;
		  double b2 = cx2 - cx1;
		  double c2 = cx1 * cy2 - cx2 * cy1;
	        
		  double D = a1*b2 - a2*b1;

		  if(D != 0)
		    {
		      double xi = (b1*c2 - b2*c1)/D; //cout<<xi<<endl;
		      double yi = (a2*c1 - a1*c2)/D; //cout<<yi<<endl;
		      hXOffsetFromMoller->Fill(xi);
		      hYOffsetFromMoller->Fill(yi);
		      hhMollerCenter->Fill(xi, yi);
		  
		      if( (symmetricTheta1 > 1.2) && (symmetricTheta2>1.2) )
			{
			  hXOffsetFromSymmetricMoller->Fill(xi);
			  hYOffsetFromSymmetricMoller->Fill(yi);
			  hhSymmetricMollerCenter->Fill(xi, yi);
			  hSymmetricMollerTheta->Fill(symmetricTheta1);
			  hSymmetricMollerTheta->Fill(symmetricTheta2);
			}
		    }

		}
	      else
		{
		  px1 = cx1 = gem1[0].x;
		  py1 = cy1 = gem1[0].y;
		  px2 = cx2 = gem2[0].x;
		  py2 = cy2 = gem2[0].y;

		}
            }  
	  }
      }
    }
  
  else if ( (gem2.size() == 2)&&(gem1.size() == 2) )
    {
      /*
       * Note:
       *     coordinates returned from class ZeroHitDecoder or OnlineHitDecoder
       *     is in chamber coordinate. 
       *     Total length in X side: 550.4mm
       *
       *     Need to find out the overlap region. 51.2mm overlapping. 
       *     checked the gerber file. only the hole area is overlapping.
       *     overlapping x region: x >  550.4/2  - 51.2 (OverlapStart)
       */
      /*
	if( (gem1[0].x>OverlapStart) && (gem2[0].x>OverlapStart)  &&  (gem1[1].x>OverlapStart) && (gem2[1].x>OverlapStart) )
	{
	if(  ( (-gem2[0].x+O_Transfer)*(-gem2[1].x+O_Transfer) < 0 ) && ( gem2[0].y*gem2[1].y < 0 ) ) 
	{
	double temp = 0;
	double slope1 = 0;
	double slope2 = 0;
	//1st electron
	theta = TMath::Sqrt( (O_Transfer-gem2[1].x)*(O_Transfer-gem2[1].x) + gem2[1].y*gem2[1].y) / z_gem1;
	theta = TMath::ATan(theta);
	slope1 = -gem2[1].y/(O_Transfer-gem2[1].x);
	slope1 = TMath::ATan(slope1);
	//hThetaDistributionMoller->Fill(theta*180.0/PI);
	temp = theta;
	//2nd electron
	theta = TMath::Sqrt( (O_Transfer-gem2[0].x)*(O_Transfer-gem2[0].x) + gem2[0].y*gem2[0].y) / z_gem2;
	theta = TMath::ATan(theta);
	slope2 = gem2[0].y/(gem2[0].x-O_Transfer);
	slope2 = TMath::ATan(slope2);
	if(theta > temp) {hThetaDistributionMollerLarge->Fill(theta); hThetaDistributionMollerSmall->Fill(temp);}
	else {hThetaDistributionMollerLarge->Fill(temp); hThetaDistributionMollerSmall->Fill(theta);}
	theta += temp;
	hThetaDistributionMoller->Fill(theta*180.0/PI);
	hLinearDeviationMoller->Fill((slope1-slope2)*180.0/PI);
	}
	}
      */
    }

}

template<class T> void GEMPhysHandler::ProcessMollerAfterCorrection(T * hit_decoder)
{
  // theta distribution
  double z_gem = 5260; //mm
  double theta = 0;
  double q_square = 0;
  double top = 0.0;
  double bottom = 0.0; //q_square = top/bottom
  double e0 = beamEnergy; //GeV
  
  vector<PRadGEMCluster> gem1, gem2;

  // GetCuster2DPositionBeamLine() returns a plane on GEM2, 
  // so the z coordinate is z_gem = z_gem2 = 5260mm
  hit_decoder->GetClusterHyCalPlusMode(gem1, gem2);

  //GeometryMollerRing(gem1, y1, x2, y2);
 
  int s = HyCalGEMPosMatch(gem1, gem2, hycal_cluster);

  if(s == 0) return ;
  //Moller events
  //theta distribution for moller

  // origin offset from beam line
  // z coordinate: GEM2 : 5260mm
  double xoffset_beam = 1.631;
  double yoffset_beam = 0.366;


  // Moller events selection
  //     require them in different quadrants, very rough
 
  if( (gem1.size() == 2)&&(gem2.size() == 0)) 
    { 
      if(  ( (gem1[0]).x*(gem1[1]).x < 0 ) && ( gem1[0].y*gem1[1].y < 0 ) ) 
	{
	  //beamline correction
	  gem1[0].x -= xoffset_beam;
	  gem1[1].x -= xoffset_beam;
	  gem1[0].y -= yoffset_beam;
	  gem1[1].y -= yoffset_beam;
	  if( (gem1[0].energy+gem1[1].energy) < totalEnergyDeposit) return;
	  double temp = 0;
	  double slope1 = 0;
	  double slope2 = 0;
	  //1st electron
	  theta = TMath::Sqrt((gem1[0].x)*(gem1[0].x) + gem1[0].y*gem1[0].y) / z_gem;
	  theta = TMath::ATan(theta);
	  hhEnergyVsAngleBeamLineCorrection->Fill(theta*180.0/PI, gem1[0].energy);
	  hhEnergyVsAngleMollerBeamLineCorrection->Fill(theta*180.0/PI, gem1[0].energy);
	  
	  slope1 = TMath::Abs( gem1[0].y/gem1[0].x );
	  slope1 = TMath::ATan(slope1); assert( slope1 < PI/2);
	  if( (gem1[0].x<0) && (gem1[0].y>0) ) slope1 = slope1+ PI/2; 
	  else if( (gem1[0].x<0) && (gem1[0].y<0) ) slope1 = slope1+PI; 
	  else if( (gem1[0].x>0) && (gem1[0].y <0)) slope1 = slope1 + 1.5*PI;
	  //hThetaDistributionMoller->Fill(theta*180.0/PI);
	  temp = theta;
	  //2nd electron
	  theta = TMath::Sqrt((gem1[1].x)*(gem1[1].x) + gem1[1].y*gem1[1].y) / z_gem;
	  theta = TMath::ATan(theta);
	  hhEnergyVsAngleBeamLineCorrection->Fill(theta*180.0/PI, gem1[1].energy);
	  hhEnergyVsAngleMollerBeamLineCorrection->Fill(theta*180.0/PI, gem1[1].energy);

	  slope2 = gem1[1].y / (gem1[1].x);
	  slope2 = TMath::ATan(slope2);
	  if( (gem1[1].x<0) && (gem1[1].y>0) ) slope1 = slope1+ PI/2; 
	  else if( (gem1[1].x<0) && (gem1[1].y<0) ) slope1 = slope1+PI; 
	  else if( (gem1[1].x>0) && (gem1[1].y<0)) slope1 = slope1 + 1.5*PI;

	  theta+=temp;
	  hThetaDistributionMollerBeamLineCorrection->Fill(theta*180.0/PI);
	  hLinearDeviationMollerAfterCorrection->Fill((slope1-slope2)*180.0/PI - 180.0);

	  {
	    if(px1_c != 0x270F)
	      {
	        px1_c = cx1_c;
		py1_c = cy1_c;
		px2_c = cx2_c;
		py2_c = cy2_c;
		cx1_c = gem1[0].x;
		cx2_c = gem1[1].x;
		cy1_c = gem1[0].y;
		cy2_c = gem1[1].y;
		double a1 = py1_c - py2_c;
		double b1 = px2_c - px1_c;
		double c1 = px1_c * py2_c - px2_c * py1_c;

		double a2 = cy1_c - cy2_c;
		double b2 = cx2_c - cx1_c;
		double c2 = cx1_c * cy2_c - cx2_c * cy1_c;
	        
		double D = a1*b2 - a2*b1;

		if(D != 0)
		  {
		    double xi = (b1*c2 - b2*c1)/D; //cout<<xi<<endl;
		    double yi = (a2*c1 - a1*c2)/D; //cout<<yi<<endl;
		    hXOffsetFromMollerAfterCorrection->Fill(xi);
		    hYOffsetFromMollerAfterCorrection->Fill(yi);
		    hhMollerCenterAfterCorrection->Fill(xi, yi);
		  }

	      }
	    else
	      {
	        px1_c = cx1_c = gem1[0].x;
		py1_c = cy1_c = gem1[0].y;
		px2_c = cx2_c = gem1[1].x;
		py2_c = cy2_c = gem1[1].y;
	      }
	  }
	} 
    }
  else if( (gem2.size() == 2)&&(gem1.size() == 0))
    { 
      if(  ( (gem2[0].x)*(gem2[1].x) < 0 ) && ( gem2[0].y*gem2[1].y < 0 ) ) 
	{
          //beamline correction
	  gem2[0].x -= xoffset_beam;
	  gem2[1].x -= xoffset_beam;
	  gem2[0].y -= yoffset_beam;
	  gem2[1].y -= yoffset_beam;

	  if( (gem2[0].energy+gem2[1].energy) < totalEnergyDeposit) return;
	  double temp = 0;
	  double slope1 = 0;
	  double slope2 = 0;
	  //1st electron
	  theta = TMath::Sqrt( (gem2[0].x)*(gem2[0].x) + gem2[0].y*gem2[0].y) / z_gem;
	  theta = TMath::ATan(theta);
	  hhEnergyVsAngleBeamLineCorrection->Fill(theta*180.0/PI, gem2[0].energy);
	  hhEnergyVsAngleMollerBeamLineCorrection->Fill(theta*180.0/PI, gem2[0].energy);

          slope1 = TMath::Abs( gem2[0].y/gem2[0].x );
	  slope1 = TMath::ATan(slope1); assert( slope1 < PI/2);
	  if( (gem2[0].x<0) && (gem2[0].y>0) ) slope1 = slope1+ PI/2; 
	  else if( (gem2[0].x<0) && (gem2[0].y<0) ) slope1 = slope1+PI; 
	  else if( (gem2[0].x>0) && (gem2[0].y <0)) slope1 = slope1 + 1.5*PI;
	  //hThetaDistributionMoller->Fill(theta*180.0/PI);
	  temp = theta;
	  //2nd electron
	  theta = TMath::Sqrt( (gem2[1].x)*(gem2[1].x) + gem2[1].y*gem2[1].y) / z_gem;
	  theta = TMath::ATan(theta);
	  hhEnergyVsAngleBeamLineCorrection->Fill(theta*180.0/PI, gem2[1].energy);
	  hhEnergyVsAngleMollerBeamLineCorrection->Fill(theta*180.0/PI, gem2[1].energy);

	  slope2 = gem2[1].y / (gem2[1].x);
	  slope2 = TMath::ATan(slope2);
	  if( (gem2[1].x<0) && (gem2[1].y>0) ) slope1 = slope1+ PI/2; 
	  else if( (gem2[1].x<0) && (gem2[1].y<0) ) slope1 = slope1+PI; 
	  else if( (gem2[1].x>0) && (gem2[1].y<0)) slope1 = slope1 + 1.5*PI;
	  theta+=temp;
	  hThetaDistributionMollerBeamLineCorrection->Fill(theta*180.0/PI);
	  hLinearDeviationMollerAfterCorrection->Fill((slope1-slope2)*180.0/PI - 180.0);

	  //x, y offset
	  //if(TMath::Abs( (slope1-slope2)*180.0/PI) < 1.0)
	  {
	    if(px1_c != 0x270F)
	      {
	        px1_c = cx1_c;
		py1_c = cy1_c;
		px2_c = cx2_c;
		py2_c = cy2_c;
		cx1_c = gem2[0].x;
		cx2_c = gem2[1].x;
		cy1_c = gem2[0].y;
		cy2_c = gem2[1].y;
		double a1 = py1_c - py2_c;
		double b1 = px2_c - px1_c;
		double c1 = px1_c * py2_c - px2_c * py1_c;

		double a2 = cy1_c - cy2_c;
		double b2 = cx2_c - cx1_c;
		double c2 = cx1_c * cy2_c - cx2_c * cy1_c;
	        
		double D = a1*b2 - a2*b1;

		if(D != 0)
		  {
		    double xi = (b1*c2 - b2*c1)/D; //cout<<xi<<endl;
		    double yi = (a2*c1 - a1*c2)/D; //cout<<yi<<endl;
		    hXOffsetFromMollerAfterCorrection->Fill(xi);
		    hYOffsetFromMollerAfterCorrection->Fill(yi);
		    hhMollerCenterAfterCorrection->Fill(xi, yi);
		  }

	      }
	    else
	      {
	        px1_c = cx1_c = gem2[0].x;
		py1_c = cy1_c = gem2[0].y;
		px2_c = cx2_c = gem2[1].x;
		py2_c = cy2_c = gem2[1].y;
	      }
	  }
	}
    }
  else if ( (gem2.size() == 1)&&(gem1.size() == 1) )
    {
      {
	if(  ( (gem1[0].x)*(gem2[0].x) < 0 ) && ( gem1[0].y* (gem2[0].y) < 0 ) ) 
	  { 
	    //cout<<gem1[0].x << " "<<gem2[0].x <<" "<<gem1[0].y << " "<<gem2[0].y<<endl;
            //beamline correction
	    gem1[0].x -= xoffset_beam;
	    gem2[0].x -= xoffset_beam;
	    gem1[0].y -= yoffset_beam;
	    gem2[0].y -= yoffset_beam;
	    //cout<<gem1[0].x << " "<<gem2[0].x <<" "<<gem1[0].y << " "<<gem2[0].y<<endl;
            //getchar();

	    if( (gem1[0].energy + gem2[0].energy) < totalEnergyDeposit) return; 
	    double temp = 0;
	    double slope1 = 0;
	    double slope2 = 0;
	    //1st electron
	    theta = TMath::Sqrt( (gem1[0].x)*(gem1[0].x) + gem1[0].y*gem1[0].y) / z_gem;
	    theta = TMath::ATan(theta);
	    hhEnergyVsAngleBeamLineCorrection->Fill(theta*180.0/PI, gem1[0].energy);
	    hhEnergyVsAngleMollerBeamLineCorrection->Fill(theta*180.0/PI, gem1[0].energy);

	    slope1 = gem1[0].y / (gem1[0].x);
	    slope1 = TMath::ATan(slope1);
	    if( (gem1[0].x<0) && (gem1[0].y>0) ) slope1 = slope1+ PI/2; 
	    else if( (gem1[0].x<0) && (gem1[0].y<0) ) slope1 = slope1+PI; 
	    else if( (gem1[0].x>0) && (gem1[0].y<0)) slope1 = slope1 + 1.5*PI;
	    //hThetaDistributionMoller->Fill(theta*180.0/PI);
	    temp = theta;
	    //2nd electron
	    theta = TMath::Sqrt( (gem2[0].x)*(gem2[0].x) + gem2[0].y*gem2[0].y) / z_gem;
	    theta = TMath::ATan(theta);
	    hhEnergyVsAngleBeamLineCorrection->Fill(theta*180.0/PI, gem2[0].energy);
	    hhEnergyVsAngleMollerBeamLineCorrection->Fill(theta*180.0/PI, gem2[0].energy);

	    slope2 = gem2[0].y / (gem2[0].x);
	    slope2 = TMath::ATan(slope2);
	    if( (gem2[0].x<0) && (gem2[0].y>0) ) slope1 = slope1+ PI/2; 
	    else if( (gem2[0].x<0) && (gem2[0].y<0) ) slope1 = slope1+PI; 
	    else if( (gem2[0].x>0) && (gem2[0].y<0)) slope1 = slope1 + 1.5*PI;
	    theta+=temp;
	    hThetaDistributionMollerBeamLineCorrection->Fill(theta*180.0/PI);
	    hLinearDeviationMollerAfterCorrection->Fill((slope1-slope2)*180.0/PI - 180.0);

	    //x, y offset
	    //if(TMath::Abs( (slope1-slope2)*180.0/PI) < 1.0)
	    {
	      if(px1_c != 0x270F)
		{
		  px1_c = cx1_c;
		  py1_c = cy1_c;
		  px2_c = cx2_c;
		  py2_c = cy2_c;
		  cx1_c = gem1[0].x;
		  cx2_c = gem2[0].x;
		  cy1_c = gem1[0].y;
		  cy2_c = gem2[0].y;
		  double a1 = py1_c - py2_c;
		  double b1 = px2_c - px1_c;
		  double c1 = px1_c * py2_c - px2_c * py1_c;

		  double a2 = cy1_c - cy2_c;
		  double b2 = cx2_c - cx1_c;
		  double c2 = cx1_c * cy2_c - cx2_c * cy1_c;
	        
		  double D = a1*b2 - a2*b1;

		  if(D != 0)
		    {
		      double xi = (b1*c2 - b2*c1)/D; 
		      double yi = (a2*c1 - a1*c2)/D;
		      //cout<<xi<<" "<<yi<<" D :"<<D<<endl;
		      //getchar();
		      hXOffsetFromMollerAfterCorrection->Fill(xi);
		      hYOffsetFromMollerAfterCorrection->Fill(yi);
		      hhMollerCenterAfterCorrection->Fill(xi, yi);
		    }

		}
	      else
		{
		  px1_c = cx1_c = gem1[0].x;
		  py1_c = cy1_c = gem1[0].y;
		  px2_c = cx2_c = gem2[0].x;
		  py2_c = cy2_c = gem2[0].y;
		}
            }  
	  }
      }
    }
}

template<class T> void GEMPhysHandler::CharactorizeGEM(T * hit_decoder)
{
  // theta distribution
  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm
  double theta = 0;
  double q_square = 0;
  double top = 0.0;
  double bottom = 0.0; //q_square = top/bottom
  double e0 = beamEnergy; //GeV

  vector<PRadGEMCluster> gem1, gem2;
  hit_decoder->GetClusterGEM(gem1, gem2);

  int nbClusterPerEvent; 
  int n1 = gem1.size();
  if(n1>0) nbClusterPerEvent = n1;
  nbClusterPerEvent += gem2.size();

  if(gem1.size() > 0)
    {
      for(int i=0;i<gem1.size();i++)
	{
	  theta = TMath::Sqrt( (gem1[i].x-O_Transfer)*(gem1[i].x-O_Transfer) + gem1[i].y*gem1[i].y) / z_gem1;
	  theta = TMath::ATan(theta); 
	  hThetaDistribution->Fill(theta*180.0/PI);

	  hhClusterDist[0] -> Fill( gem1[i].x, gem1[i].y);
	  hhClusterDist_pRad[0]->Fill( gem1[i].x, (-1.0)*gem1[i].y );
	  hhChargeRatio[0] -> Fill(gem1[i].x_charge, gem1[i].y_charge);
	  hChargeRatio[0] -> Fill( gem1[i].x_charge/gem1[i].y_charge);
	  hClusterADCDistXSide[0]->Fill( gem1[i].x_charge);
	  hClusterADCDistYSide[0]->Fill( gem1[i].y_charge);
	  if( gem1[i].x_charge<0)cout<<"x1 negative charge???"<<endl;;
	  if( gem1[i].y_charge<0)cout<<"y1 negative charge???"<<endl;;
	  if( gem1[i].x_charge<100 || gem1[i].y_charge<100 )
	    {
	      //cout<<"Warning!!! unusual cluster adcs..."<<endl;
	      //cout<<"GEM1: charge: "<<charge_gem1[i].x<<"  "<<charge_gem1[i].y<<endl;
	    }
	}
    }

  if(gem2.size() > 0)
    {
      for(int i=0;i<gem2.size();i++)
	{
	  //nb of clusters per event
	  if(gem2[i].x >= OverlapStart) 
	    {
	      theta = TMath::Sqrt( (O_Transfer-gem2[i].x)*(O_Transfer-gem2[i].x) + gem2[i].y*gem2[i].y ) / z_gem2;
	      theta = TMath::ATan(theta);
	      hThetaDistribution->Fill(theta*180.0/PI);
	      nbClusterPerEvent--;
	    }

	  hhClusterDist[1] -> Fill( gem2[i].x, gem2[i].y);
	  hhClusterDist_pRad[1]->Fill((-1.0)* gem2[i].x, gem2[i].y );
	  hhChargeRatio[1] -> Fill(gem2[i].x_charge, gem2[i].y_charge);
	  hChargeRatio[1] -> Fill( gem2[i].x_charge/gem2[i].y_charge);
	  hClusterADCDistXSide[1]->Fill( gem2[i].x_charge);
	  hClusterADCDistYSide[1]->Fill( gem2[i].y_charge);
	  if( gem2[i].x_charge<0)cout<<"x2 negative charge???"<<endl;
	  if( gem2[i].y_charge<0)cout<<"y2 negative charge???"<<endl;
	  if( gem2[i].x_charge<100 || gem2[i].y_charge<100 )
	    {
	      //cout<<"Warning!!! unusual cluster adcs..."<<endl;
	      //cout<<"GEM2: charge: "<<charge_gem2[i].x<<"  "<<charge_gem2[i].y<<endl;
	    }
	}
    }
  hGEMClusterMul->Fill(nbClusterPerEvent);
}

void GEMPhysHandler::SavePhysResults()
{
  const char *str = config.phys_results_path.c_str();
  TFile *f = new TFile(str, "recreate");
  
  // gem charatorization
  TDirectory * gem = f->mkdir("GEM_Charactorization");
  gem->cd();
  hhClusterDist[0]->Write();
  hhClusterDist_pRad[0]->Write();
  hhChargeRatio[0]->Write();
  hChargeRatio[0]->Write();
  hClusterADCDistXSide[0]->Write();
  hClusterADCDistYSide[0]->Write();
  hNbClusterPerPlaneX[0]->Write();
  hClusterDistX[0]->Write();
  hNbClusterPerPlaneY[0]->Write();
  hClusterDistY[0]->Write();

  hhClusterDist[1]->Write();
  hhClusterDist_pRad[1]->Write();

  hhChargeRatio[1]->Write();
  hChargeRatio[1]->Write();
  hClusterADCDistXSide[1]->Write();
  hClusterADCDistYSide[1]->Write();
  hNbClusterPerPlaneX[1]->Write();
  hClusterDistX[1]->Write();
  hNbClusterPerPlaneY[1]->Write();
  hClusterDistY[1]->Write();

  // gem in combination with hycal
  TDirectory * gem_hycal = f->mkdir("GEM_HyCal_charatorization");
  gem_hycal->cd();
  hhGEMClusterMap->Write();
  hhGEMClusterMapSingleCluster->Write();
  hhGEMClusterMapMaxTwoCluster->Write();
  hhHyCalClusterMap->Write();
  hHyCalEnergy->Write();
  hHyCalEnergyEp->Write();
  hHyCalEnergyMoller->Write();

  hGEMClusterMul->Write();
  hHyCalClusterMul->Write();

  // theta distribution
  TDirectory * theta = f->mkdir("Theta_Distribution_coplanarity");
  theta->cd();
  hThetaDistribution->Write();
  hThetaDistributionEp->Write();
  hThetaDistributionMoller->Write();
  hThetaDistributionMollerLarge->Write();
  hThetaDistributionMollerSmall->Write();
  hThetaDistributionMollerSlice->Write();
  hLinearDeviationMoller->Write();
  hLinearDeviationMollerAfterCorrection->Write();
  hThetaDistributionMollerBeamLineCorrection->Write();


  //q square distribution
  TDirectory * q_square = f->mkdir("Q_Square_Distribution");
  q_square->cd();
  hQSquareEp->Write();
  hQSquareEp1->Write();
  hQSquareEp2->Write();
  hQSquareEp3->Write();
  hQSquareEp4->Write();
  hQSquareEp5->Write();
  hhQSquareScattAngleEp->Write();

  // every vs angle
  TDirectory * energyangle = f->mkdir("Energy_Angle");
  energyangle->cd();
  hhEnergyVsAngle->Write();
  hhEnergyVsAngleMoller->Write();
  hhEnergyVsAngleEp->Write();
  hhEnergyVsAngleBeamLineCorrection->Write();
  hhEnergyVsAngleMollerBeamLineCorrection->Write();
  hhEnergyVsAngleHyCal2CGEM1C->Write();
  
  // hycal two cluster gem only 1
  TDirectory * hycal2gem1 = f->mkdir("HyCal2Cluster_GEM1Cluster");
  hycal2gem1->cd();
  hhGEMClusterMapHyCal2GEM1->Write();
  hhHyCalClusterMapHyCal2GEM1->Write();
  hhHyCalClusterMap2ClusterBeforeMatch->Write();
  hhGEMClusterMap2ClusterAfterMatch->Write();

  // offset
  TDirectory * offset = f -> mkdir("Offsets");
  offset->cd();
  hXOffsetFromMoller->Write();
  hYOffsetFromMoller->Write();
  hXOffsetFromMollerAfterCorrection->Write();
  hYOffsetFromMollerAfterCorrection->Write();
  hXOffsetFromSymmetricMoller->Write();
  hYOffsetFromSymmetricMoller->Write();
  hXOffsetGEM->Write();
  hYOffsetGEM->Write();
  hXOffsetGEMAfterCorrection->Write();
  hYOffsetGEMAfterCorrection->Write();


  // moller ring
  TDirectory *moller_ring = f->mkdir("Moller_Ring");
  moller_ring->cd();
  hhMollerCenter->Write();
  hhMoller2DRing->Write();
  hhMoller2DRing2->Write();
  hhMoller2DRingBeforeMatch->Write();
  hhMoller2DRingBeforeMatch2->Write();
  hhMollerCenterAfterCorrection->Write();
  hhSymmetricMollerCenter->Write();
  hSymmetricMollerTheta->Write();
  
  // hycal gem match
  TDirectory * hycal_match = f->mkdir("Hycal_GEM_Match");
  hycal_match->cd();
  hXDiffMatch->Write();
  hYDiffMatch->Write();
  hhXDiffMatchVsEnergy->Write();
  hhYDiffMatchVsEnergy->Write();
  hNbPointsMatch->Write();

  // un-relative
  TDirectory * u = f->mkdir("time_correlation");
  u->cd();
  hhTimeCorrelation->Write();
  hTimeDiff->Write();

  f->Write();
  f->Save();
  //f->Close();
}

int GEMPhysHandler::GEMHyCalPosMatch(int ngem, vector<PRadGEMCluster> &gem, vector<HyCalHit> *pHHit)
{
  if(gem.size() == 0)
    {
      pHHit->clear();
      return 0;
    }
  if(pHHit->size() == 0) 
    {
      gem.clear();
      return 0;
    }

  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm
  double z_gem = 0;
  double z_hycal = 5820; //mm
  double res = 10; // a larger range, 60mm

  if( ngem==0 ) 
    z_gem = z_gem1; 
  else if(ngem == 1) 
    z_gem = z_gem2;
  else
    {
      cout<<"GEMIndex: 0 or 1..."<<endl;
      return 0;
    }
  
  vector<PRadGEMCluster> res_gem;

  int n = gem.size();
  int nh = pHHit->size();
  for(int i=0;i<n;i++)
    {
      double r_gem = TMath::Sqrt( gem[i].x*gem[i].x + gem[i].y*gem[i].y );
      double r_hycal = r_gem/z_gem * z_hycal;
      double x_hycal = r_hycal*gem[i].x/r_gem;
      double y_hycal = r_hycal*gem[i].y/r_gem;

      for(int j=0;j<nh;j++)
	{
	  if( TMath::Abs((pHHit->at(j).x)-x_hycal) < res)
	    {
	      if( TMath::Abs((pHHit->at(j).y)-y_hycal) < res)
		{
		  gem[i].energy = pHHit->at(j).E;
		  res_gem.push_back(gem[i]);
		}
	    }
	}
    }

  if(res_gem.size() > 0)
    {
      gem.clear();
      gem = res_gem;
    
      return res_gem.size();
    }
  else
    {
      gem.clear();
      return 0;

    }
}

int GEMPhysHandler::HyCalGEMPosMatch( vector<PRadGEMCluster> &gem1, vector<GEMClusterStruct> &gem2, vector<HyCalHit> *pHHit)
{
  int nhits_gem1 = gem1.size();
  int nhits_gem2 = gem2.size();
  int nhits_hycal = pHHit->size();

  if( (gem1.size() == 0) && (gem2.size() == 0))
    {
 
      pHHit->clear();
      return 0;
    }
  if(pHHit->size() == 0) 
    {
      gem1.clear();
      gem2.clear();
 
      return 0;
    }

  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm
  double z_gem = 0;
  double z_hycal = 5820; //mm
  double res = 10; // a larger range, 60mm

  vector<PRadGEMCluster> res_gem1;
  vector<PRadGEMCluster> res_gem2;

  int nh = pHHit->size();
  for(int i=0;i<nh;i++)
    {
      int dx = res;
      int dy = res;

      int match_gem1 = 0;
      int match_gem2 = 0;
      float m_index = 0;
      float m_e = 0;

      // search GEM1
      int n = gem1.size();
      if(n > 0)
	{
	  double x_hycal = (pHHit->at(i).x) *z_gem1/z_hycal;
	  double y_hycal = (pHHit->at(i).y) *z_gem1/z_hycal;
	  for(int j=0;j<n;j++)
	    {
	      if( (TMath::Abs(gem1[j].x-x_hycal) < dx) && ( TMath::Abs(gem1[j].y-y_hycal) < dy) )
		{
		  dx = TMath::Abs(gem1[j].x-x_hycal);
		  dy = TMath::Abs(gem1[j].y-y_hycal); 
		  m_index = j;
		  m_e = pHHit->at(i).E;
		  match_gem1 = 1;
		  match_gem2 = 0;
		  hXDiffMatch->Fill(gem1[j].x - x_hycal);
		  hYDiffMatch->Fill(gem1[j].y - y_hycal);
		  hhXDiffMatchVsEnergy->Fill(m_e, gem1[j].x - x_hycal);
		  hhYDiffMatchVsEnergy->Fill(m_e, gem1[j].y - y_hycal);
		}
	    }
	}

      //continue search GEM2
      n = gem2.size();
      if(n>0)
	{
	  double x_hycal = (pHHit->at(i).x) *z_gem2/z_hycal;
	  double y_hycal = (pHHit->at(i).y) *z_gem2/z_hycal;
	  for(int j=0;j<n;j++)
	    {
	      if( (TMath::Abs(gem2[j].x-x_hycal) < dx) && ( TMath::Abs(gem2[j].y-y_hycal) < dy) )
		{
		  dx = TMath::Abs(gem2[j].x-x_hycal);
		  dy = TMath::Abs(gem2[j].x-y_hycal); 
		  m_index = j;
		  m_e = pHHit->at(i).E;
		  match_gem2 = 1;
		  match_gem1 = 0;
		  hXDiffMatch->Fill(gem2[j].x - x_hycal);
		  hYDiffMatch->Fill(gem2[j].y - y_hycal);
		  hhXDiffMatchVsEnergy->Fill(m_e, gem2[j].x - x_hycal);
		  hhYDiffMatchVsEnergy->Fill(m_e, gem2[j].y - y_hycal);
		}
	    }
	}
      //after searching
      if( (match_gem1 == 1) && (match_gem2 == 1) ) cout<<"GEMPhysHandler::HyCalGEMPosMatch(): error...."<<endl;
      if( match_gem1 == 1 ) { gem1[m_index].energy = m_e; res_gem1.push_back(gem1[m_index]);}
      if( match_gem2 == 1 ) { gem2[m_index].energy = m_e; res_gem2.push_back(gem2[m_index]);}
    }
  
  gem1.clear();
  gem2.clear();
  if(res_gem1.size() > 0)
    {
      gem1 = res_gem1;
    
    }
  if(res_gem2.size() > 0)
    {
      gem2 = res_gem2;
    }

  if(nhits_hycal == 2)
    {
      for(int i=0;i<nhits_hycal;i++)
	{
	  hhHyCalClusterMap2ClusterBeforeMatch->Fill(pHHit->at(i).x, pHHit->at(i).y);
	}
    }

  if( gem1.size() + gem2.size() == 2 )
    {
      // to detect gem moller efficiency
      // see after match , how many two cluster events in gem
      for(int i = 0;i<res_gem1.size();i++)
	{
	  hhGEMClusterMap2ClusterAfterMatch->Fill(res_gem1[i].x, res_gem1[i].y);
	}
      for(int i = 0;i<res_gem2.size();i++)
	{
	  hhGEMClusterMap2ClusterAfterMatch->Fill(res_gem2[i].x, res_gem2[i].y);
	}
    }

  if( (gem1.size() + gem2.size()) > 0)
    {
      if( (nhits_gem1 == 1) && (nhits_gem2==0) &&  (nhits_hycal >=2) )
	{
	  double theta = TMath::Sqrt( (res_gem1[0].x)*(res_gem1[0].x)+(res_gem1[0].y)*(res_gem1[0].y)  )/z_gem1 *180 /PI;
	  hhEnergyVsAngleHyCal2CGEM1C->Fill(theta, gem1[0].energy);
	  hhGEMClusterMapHyCal2GEM1->Fill( res_gem1[0].x, res_gem1[0].y);
	  for(int i=0;i<nhits_hycal;i++)
	    {
	      hhHyCalClusterMapHyCal2GEM1->Fill(pHHit->at(i).x, pHHit->at(i).y);
	    }
	}
      if( (nhits_gem2 == 1) && (nhits_gem1==0) &&  (nhits_hycal >=2) )
	{
	  double theta = TMath::Sqrt( (res_gem2[0].x)*(res_gem2[0].x)+(res_gem2[0].y)*(res_gem2[0].y)  )/z_gem2 *180 /PI;
	  hhEnergyVsAngleHyCal2CGEM1C->Fill(theta, gem2[0].energy);
	  hhGEMClusterMapHyCal2GEM1->Fill( res_gem2[0].x, res_gem2[0].y);
	  for(int i=0;i<nhits_hycal;i++)
	    {
	      hhHyCalClusterMapHyCal2GEM1->Fill(pHHit->at(i).x, pHHit->at(i).y);
	    }

	}
    }

  hNbPointsMatch->Fill( res_gem1.size() + res_gem2.size() );

#ifndef CLEANUP
  return res_gem1.size() + res_gem2.size();
#endif

#ifdef CLEANUP
  if ( res_gem1.size() + res_gem2.size() > 2)
    return 0;
  else 
    return res_gem1.size() + res_gem2.size();
#endif
}

void GEMPhysHandler::GeometryMollerRing(vector<PRadGEMCluster> &gem1, vector<GEMClusterStruct>& gem2)
{
  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm
  double theta = 0;

  // moller 2d ring
  float ThetaSmall = 0.7/180.0 * PI;
  float ThetaLarge = 0.8/180.0 * PI;
  float ThetaSmall2 = 1.0/180.0 * PI;
  float ThetaLarge2 = 1.1/180.0 * PI;

  // Moller events selection
  //     require them in different quadrants, very rough
 
  if( (gem1.size() == 2)&&(gem2.size() == 0)) 
    { 
      if(  ( (gem1[0].x)*(gem1[1].x) < 0 ) && ( gem1[0].y*gem1[1].y < 0 ) ) 
	{
	  double temp = 0;
	  //1st electron
	  theta = TMath::Sqrt((gem1[0].x)*(gem1[0].x) + gem1[0].y*gem1[0].y) / z_gem1;
	  theta = TMath::ATan(theta);
	  
	  temp = theta;
	  //2nd electron
	  theta = TMath::Sqrt((gem1[1].x)*(gem1[1].x) + gem1[1].y*gem1[1].y) / z_gem1;
	  theta = TMath::ATan(theta);

	  //2d ring
	  if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge)   ) )
	    {
	      hhMoller2DRingBeforeMatch->Fill(gem1[0].x, gem1[0].y);
	      hhMoller2DRingBeforeMatch->Fill(gem1[1].x, gem1[1].y);
	    }
	  if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2)   ) )
	    {
	      hhMoller2DRingBeforeMatch2->Fill(gem1[0].x, gem1[0].y);
	      hhMoller2DRingBeforeMatch2->Fill(gem1[1].x, gem1[1].y);
	    }
	} 
    }
  else if( (gem2.size() == 2)&&(gem1.size() == 0))
    { 
      if(  ( (gem2[0].x)*(gem2[1].x) < 0 ) && ( gem2[0].y*gem2[1].y < 0 ) ) 
	{
	  double temp = 0;
	  //1st electron
	  theta = TMath::Sqrt( (gem2[0].x)*(gem2[0].x) + gem2[0].y*gem2[0].y) / z_gem2;
	  theta = TMath::ATan(theta);

	  temp = theta;
	  //2nd electron
	  theta = TMath::Sqrt( (gem2[1].x)*(gem2[1].x) + gem2[1].y*gem2[1].y) / z_gem2;
	  theta = TMath::ATan(theta);

	  //2d ring
	  if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge)   ) )
	    {
	      hhMoller2DRingBeforeMatch->Fill(gem2[0].x, gem2[0].y);
	      hhMoller2DRingBeforeMatch->Fill(gem2[1].x, gem2[1].y);
	    }
	  if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2)   ) )
	    {
	      hhMoller2DRingBeforeMatch2->Fill(gem2[0].x, gem2[0].y);
	      hhMoller2DRingBeforeMatch2->Fill(gem2[1].x, gem2[1].y);
	    }
	}
    }
  else if ( (gem2.size() == 1)&&(gem1.size() == 1) )
    {
      {
	if(  ( (gem1[0].x)*(gem2[0].x) < 0 ) && ( gem1[0].y* (gem2[0].y) < 0 ) ) 
	  { 
	    double temp = 0;
	    //1st electron
	    theta = TMath::Sqrt( (gem1[0].x)*(gem1[0].x) + gem1[0].y*gem1[0].y) / z_gem1;
	    theta = TMath::ATan(theta);

	    temp = theta;
	    //2nd electron
	    theta = TMath::Sqrt( (gem2[0].x)*(gem2[0].x) + gem2[0].y*gem2[0].y) / z_gem2;
	    theta = TMath::ATan(theta);

	    //2d ring
	    if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge)   ))
              {
	        hhMoller2DRingBeforeMatch->Fill(gem2[0].x, gem2[0].y);
	        hhMoller2DRingBeforeMatch->Fill(gem1[0].x, gem1[0].y);
	      }
	    if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2)   ))
              {
	        hhMoller2DRingBeforeMatch2->Fill(gem2[0].x, gem2[0].y);
	        hhMoller2DRingBeforeMatch2->Fill(gem1[0].x, gem1[0].y);
	      }
	  }
      }
    } 
}


