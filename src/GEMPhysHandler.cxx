#include "GEMPhysHandler.h"

#include "PRadDataHandler.h"
#include "PRadReconstructor.h"
#include "PRadTDCGroup.h"
#include "PRadBenchMark.h"
#include "PRadDSTParser.h"

using namespace std;

#define PI 3.1415926

GEMPhysHandler::GEMPhysHandler()
{
  pHandler = new PRadDataHandler();
  pDSTParser = new PRadDSTParser(pHandler);

  pHandler -> ReadConfig("config.txt");

  pGEMReconstruct = new PRadGEMReconstructor(pHandler);

  //totalEnergyDeposit = 1800; //MeV
  //beamEnergy = 2.147;//GeV
  totalEnergyDeposit = 100; //MeV
  beamEnergy = 1.1;//GeV

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
  hhTimeCorrelation = 
    new TH2F("hhTimeCorrelation", "HyCal vs Scin", 1000, 0, 10000, 1000, 0, 10000);
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
      //ProcessAllEvents(-1);
    }
 
  //temporary solution 
  filename = "./prad_001287.dst";
  ProcessAllEvents( -1 );

  // Save Histos
  SavePhysResults();
}

int GEMPhysHandler::ProcessAllEvents(int evtID )
{
  cout<<"Process File:  "<<filename<<endl;
  int entry = 0;

  PRadBenchMark timer;
  
  pDSTParser->OpenInput(filename.c_str());

  while(pDSTParser->Read() && entry < 300000)
    {
      if( pDSTParser->EventType() == PRad_DST_Event) 
	{
	  auto event = pDSTParser -> GetEvent();
          auto Gem_cluster = pGEMReconstruct->CoarseGEMReconstruct(event);
          auto Hycal_cluster = pHandler->GetHyCalCluster(event);
          gem_cluster = &Gem_cluster;
          hycal_cluster = &Hycal_cluster;

	  ++entry;
	  Analyzer();
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

  gem1 = pGEMReconstruct->GEMClusteringLocal(0);
  gem2 = pGEMReconstruct->GEMClusteringLocal(1);

  hgem1 = pGEMReconstruct->GEMClusteringHyCal(0);
  hgem2 = pGEMReconstruct->GEMClusteringHyCal(1);

  //online_hit.FillHistos(hNbClusterPerPlaneX, hNbClusterPerPlaneY, hClusterDistX, hClusterDistY);

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

  CharactorizeGEM();
  ProcessEp();
  ProcessMoller();

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

void GEMPhysHandler::ProcessEp()
{
  // theta distribution
  double theta = 0;
  double q_square = 0;
  double top = 0.0;
  double bottom = 0.0; //q_square = top/bottom
  double e0 = beamEnergy; //GeV

  vector<PRadGEMCluster> gem;
  gem  = *gem_cluster;
  
  //ep events
  //theta distribution for single cluster
  if( gem.size() == 1 ) 
    { 
      // Ep energy cut
      if(gem[0].energy < totalEnergyDeposit ) return;
      double z_gem = gem[0].z;
      theta = TMath::Sqrt((gem[0].x)*(gem[0].x) + gem[0].y*gem[0].y) / z_gem;
      theta = TMath::ATan(theta);
      double theta_d = theta*180.0/PI;
      hThetaDistributionEp->Fill(theta_d);
      hhEnergyVsAngle->Fill(theta_d, gem[0].energy);
      hhEnergyVsAngleEp->Fill(theta_d, gem[0].energy);

      //q suqare
      top = 4.0*beamEnergy*beamEnergy*TMath::Sin(theta/2)*TMath::Sin(theta/2);
      bottom = 1+ (2*beamEnergy/0.938)*TMath::Sin(theta/2)*TMath::Sin(theta/2);
      q_square = top/bottom;
      hQSquareEp->Fill(q_square);
      hhQSquareScattAngleEp -> Fill(theta_d, q_square);

      if( (theta_d > 0.7) && (theta_d<0.8) ) 
	hQSquareEp1->Fill(q_square);

      if( (theta_d > 1.0) && (theta_d<1.1) ) 
	hQSquareEp2->Fill(q_square);

      if( (theta_d > 1.5) && (theta_d<1.6) ) 
	hQSquareEp3->Fill(q_square);

      if( (theta_d > 2.0) && (theta_d<2.1) ) 
	hQSquareEp4->Fill(q_square);

      if(theta_d > 2.2) 
	hQSquareEp5->Fill(q_square);
    }
}

void GEMPhysHandler::ProcessMoller()
{
  // theta distribution
  double theta = 0;
  double q_square = 0;
  double top = 0.0;
  double bottom = 0.0; //q_square = top/bottom
  double e0 = beamEnergy; //GeV
  
  vector<PRadGEMCluster> gem1, gem2;
  gem1 = pGEMReconstruct->GEMClusteringHyCal(0);
  gem2 = pGEMReconstruct->GEMClusteringHyCal(0);

  GeometryMollerRing(gem1, gem2);

  vector<PRadGEMCluster> gem;
  gem = *gem_cluster;
   
  // moller 2d ring
  float ThetaSmall = 0.7/180.0 * PI;
  float ThetaLarge = 0.8/180.0 * PI;
  float ThetaSmall2 = 1.0/180.0 * PI;
  float ThetaLarge2 = 1.1/180.0 * PI;

  // Moller events selection
  if( gem.size() == 2) 
    { 
      if( (gem[0].energy+gem[1].energy) < totalEnergyDeposit) return;
      double temp = 0;
      double slope1 = 0;
      double slope2 = 0;
      //1st electron
      theta = TMath::Sqrt((gem[0].x)*(gem[0].x) + gem[0].y*gem[0].y) / gem[0].z;
      theta = TMath::ATan(theta);
      hhEnergyVsAngle->Fill(theta*180.0/PI, gem[0].energy);
      hhEnergyVsAngleMoller->Fill(theta*180.0/PI, gem[0].energy);
	  
      slope1 = TMath::Abs( gem[0].y/gem[0].x );
      slope1 = TMath::ATan(slope1); assert( slope1 < PI/2);
      if( (gem[0].x<0) && (gem[0].y>0) ) slope1 = slope1+ PI/2; 
      else if( (gem[0].x<0) && (gem[0].y<0) ) slope1 = slope1+PI; 
      else if( (gem[0].x>0) && (gem[0].y <0)) slope1 = slope1 + 1.5*PI;
      //hThetaDistributionMoller->Fill(theta*180.0/PI);
      temp = theta;
      //2nd electron
      theta = TMath::Sqrt((gem[1].x)*(gem[1].x) + gem[1].y*gem[1].y) / gem[1].z;
      theta = TMath::ATan(theta);
      hhEnergyVsAngle->Fill(theta*180.0/PI, gem[1].energy);
      hhEnergyVsAngleMoller->Fill(theta*180.0/PI, gem[1].energy);

      slope2 = gem[1].y / (gem[1].x);
      slope2 = TMath::ATan(slope2);
      if( (gem[1].x<0) && (gem[1].y>0) ) slope1 = slope1+ PI/2; 
      else if( (gem[1].x<0) && (gem[1].y<0) ) slope1 = slope1+PI; 
      else if( (gem[1].x>0) && (gem[1].y <0)) slope1 = slope1 + 1.5*PI;
      if(theta > temp) 
	{
	  hThetaDistributionMollerLarge->Fill(theta*180.0/PI); 
	  hThetaDistributionMollerSmall->Fill(temp*180.0/PI);
	  if( ((temp*180.0/PI)>0.52) && ((temp*180.0/PI)<0.62) ) 
	    hThetaDistributionMollerSlice->Fill( theta*180.0/PI);
	}
      else 
	{
	  hThetaDistributionMollerLarge->Fill(temp*180.0/PI); 
	  hThetaDistributionMollerSmall->Fill(theta*180.0/PI);
	  if( ((theta*180.0/PI)>0.54) && ((theta*180.0/PI)<0.56) ) 
	    hThetaDistributionMollerSlice->Fill( temp*180.0/PI);
	}

      //2d ring
      if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge) ) )
	{
	  hhMoller2DRing->Fill(gem[0].x, gem[0].y);
	  hhMoller2DRing->Fill(gem[1].x, gem[1].y);
	}
      if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2) ) )
	{
	  hhMoller2DRing2->Fill(gem[0].x, gem[0].y);
	  hhMoller2DRing2->Fill(gem[1].x, gem[1].y);
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
	    cx1 = gem[0].x;
	    cx2 = gem[1].y;
	    cy1 = gem[0].y;
	    cy2 = gem[1].y;

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
	    px1 = cx1 = gem[0].x;
	    py1 = cy1 = gem[0].y;
	    px2 = cx2 = gem[1].x;
	    py2 = cy2 = gem[1].y;

	  }
      } 
    }
}

void GEMPhysHandler::CharactorizeGEM()
{
  // theta distribution
  double theta = 0;
  double q_square = 0;
  double top = 0.0;
  double bottom = 0.0; //q_square = top/bottom
  double e0 = beamEnergy; //GeV

  vector<PRadGEMCluster> gem1, gem2;
  gem1 = pGEMReconstruct->GEMClusteringLocal(0);
  gem2 = pGEMReconstruct->GEMClusteringLocal(1);

  int nbClusterPerEvent; 
  int n1 = gem1.size();
  if(n1>0) nbClusterPerEvent = n1;
  nbClusterPerEvent += gem2.size();

  if(gem1.size() > 0)
    { 
      for(int i=0;i<gem1.size();i++)
	{
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

void GEMPhysHandler::GeometryMollerRing(vector<PRadGEMCluster> &gem1, vector<PRadGEMCluster>& gem2)
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
	  if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge) ) )
	    {
	      hhMoller2DRingBeforeMatch->Fill(gem1[0].x, gem1[0].y);
	      hhMoller2DRingBeforeMatch->Fill(gem1[1].x, gem1[1].y);
	    }
	  if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2) ) )
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
	  if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge) ) )
	    {
	      hhMoller2DRingBeforeMatch->Fill(gem2[0].x, gem2[0].y);
	      hhMoller2DRingBeforeMatch->Fill(gem2[1].x, gem2[1].y);
	    }
	  if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2) ) )
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
	    if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge)))
              {
	        hhMoller2DRingBeforeMatch->Fill(gem2[0].x, gem2[0].y);
	        hhMoller2DRingBeforeMatch->Fill(gem1[0].x, gem1[0].y);
	      }
	    if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2)))
              {
	        hhMoller2DRingBeforeMatch2->Fill(gem2[0].x, gem2[0].y);
	        hhMoller2DRingBeforeMatch2->Fill(gem1[0].x, gem1[0].y);
	      }
	  }
      }
    } 
}
