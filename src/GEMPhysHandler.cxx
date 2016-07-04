#include "GEMPhysHandler.h"

#include "GEMOnlineHitDecoder.h"
#include "GEMZeroHitDecoder.h"

#include "PRadDataHandler.h"
#include "PRadReconstructor.h"
#include "PRadTDCGroup.h"
#include "PRadBenchMark.h"

//#include "PRadReconstructor.h"
using namespace std;
using namespace evio;

// GEM Efficiency Macros
#define TDC_CHECK_EFFICIENCY
//#undef  TDC_CHECK_EFFICIENCY

#define ZEROSUPPRESSION
#undef ZEROSUPPRESSION

//#define HyCalTimingCut
//#undef HyCalTimingCut

#define CLEANUP
//#undef CLEANUP

#define PI 3.1415926

GEMPhysHandler::GEMPhysHandler()
{
  fRawDecoder = 0;

  nElectron = 0;
  nElectron_126=0;
  nElectron_127=0;

  nScinEvents = 0;
  nHyCalEvents = 0;
  nTotalEvents = 0;  // to get how many events in total, this to compute photon conversion probability
  neff = 0;

  mAPVRawHistos.clear();
  mAPVRawTSs.clear();
  vSRSSingleEventData.clear();
  vSRSZeroEventData.clear();
  FECs.clear();

  PRDMapping *mapping = PRDMapping::GetInstance();
  FECs = mapping->GetBankIDSet();
  
  string pedestal_file = config.GetLoadPedPath();
  ped = new GEMPedestal(pedestal_file);
  ped -> LoadPedestal( );
  //for debug
  cout<<"FECs Found:  ";
  set<int>::iterator it;
  for(it=FECs.begin(); it!=FECs.end(); ++it)
    {
      cout<<(*it)<<"  ";
    }
  cout<<endl;

  //chao
  pHandler = new PRadDataHandler();
  pHandler->ReadTDCList("/home/xbai/w/pRad/source/PRadDecoder/config/tdc_group_list.txt");
  pHandler->ReadChannelList("/home/xbai/w/pRad/source/PRadDecoder/config/module_list.txt");
  pHandler->BuildChannelMap();
  pHandler->ReadPedestalFile("/home/xbai/w/pRad/source/PRadDecoder/config/pedestal.dat");
  pHandler->ReadCalibrationFile("/home/xbai/w/pRad/source/PRadDecoder/config/calibration.txt");
  //pHandler->EnableReconstruction();

  reconstruct = new PRadReconstructor();
  reconstruct->SetHandler(pHandler);

  totalEnergyDeposit = 1800; //MeV
  beamEnergy = 2.147;//GeV
  //compute intersection points
  px1 = 0x270F;py1=0x270F;px2=0x270F;py2=0x270F;
  cx1 = 0x270F;cy1=0x270F;cx2=0x270F;cy2=0x270F;

  px1_c = 0x270F;py1_c=0x270F;px2_c=0x270F;py2_c=0x270F;
  cx1_c = 0x270F;cy1_c=0x270F;cx2_c=0x270F;cy2_c=0x270F;
 
  //test
  hhTimeCorrelation = new TH2F("hhTimeCorrelation", "HyCal vs Scin", 1000, 0, 10000, 1000, 0, 10000);
  hTimeDiff = new TH1F("hTimeDiff", "Scin - HyCal", 1000, -5000, 5000);

  BookHistos();

  // origin transfer
  O_Transfer = 253.2;
  OverlapStart = 231.2;

}

GEMPhysHandler::~GEMPhysHandler()
{
  mAPVRawHistos.clear();
  mAPVRawTSs.clear();
  ped->Delete();
}

void GEMPhysHandler::ProcessAllFiles()
{
  //GEMConfigure config;
  int nFile = config.nFile;
  cout<<"GEMPhysHandler::ProcessAllFiles: # of Files to analyze: "<<nFile<<endl;
  for(int i=0;i<nFile;i++)
    {
      cout<<config.fileList[i]<<endl;
    }
  cout<<"TDC Cut: "<<config.TDC_Channel<<", Start TIME:  "<<config.TDC_Start<<", END TIME:  "<<config.TDC_End<<", HyCal Energy Cut:  "<<config.Hycal_Energy<<endl;

  // to compute photon conversion rate
  nTotalEvents = 0;
  neff = 0;
  nElectron = 0;


  if(config.fileList[0].find("evio.0")!= string::npos)
    {
      pHandler->InitializeByData(config.fileList[0].c_str());
    };

  for(int i=1;i<nFile;i++)
    {
      filename = config.fileList[i];
      ProcessAllEvents(-1);
      //ProcessAllEvents(40000);
    }

  // Save Histos
  SavePhysResults();

#ifdef TDC_CHECK_EFFICIENCY
  cout<<endl<<endl;
  double eff_real = (double)neff/nElectron;
  cout<<"GEM Hit: "<<neff<<endl;
  cout<<"TDC:  "<<nElectron<<endl;
  cout<<"Total Number of Events in "<<nFile<<" Files: "<<nTotalEvents<<endl;
  cout<<"------------------------------------"<<endl;
  cout<<"Overall Detector Efficiency: "<<eff_real<<endl;
  cout<<"------------------------------------"<<endl;
#endif

}

int GEMPhysHandler::ProcessAllEvents(int evtID )
{
  cout<<"Process File:  "<<filename<<endl;
  int entry = 0;
  double ntrigger = 0.0; // for detector efficiency

  PRadBenchMark timer;
  try{
    evioFileChannel chan(filename.c_str(), "r");
    chan.open();
    while(chan.read())
      {
        ntrigger +=1.0;
	nTotalEvents += 1.0;
        pHandler->Decode(chan.getBuffer());
	pHyCalHit =& reconstruct->CoarseHyCalReconstruct(pHandler->GetEventCount() - 1); 
        //pHyCalHit = &HyCalHit;
        //cout<<pHyCalHit->size()<<endl;
        hHyCalClusterMul->Fill(pHyCalHit->size());
       	for(int i=0;i<pHyCalHit->size(); i++)
	  {
	    hhHyCalClusterMap->Fill(pHyCalHit->at(i).x, pHyCalHit->at(i).y);

	    hHyCalEnergy->Fill(pHyCalHit->at(i).E);
	    if(pHyCalHit->size() == 1) hHyCalEnergyEp->Fill(pHyCalHit->at(i).E);
	    if(pHyCalHit->size() == 2) {hHyCalEnergyMoller->Fill(pHyCalHit->at(i).E);}
	  }
       
	vSRSSingleEventData.clear();
        vSRSZeroEventData.clear();

	evioDOMTree event(chan);
	evioDOMNodeListP fecEventList = event.getNodeList( isLeaf() );
	//cout<<"total number of all banks: "<<fecEventList->size()<<endl;

	evioDOMNodeList::iterator iter;

        int convert = 1;
        int HyCal_pos_match = 1;
	int HyCalTimingCut = 1;

#ifdef TDC_CHECK_EFFICIENCY
        // photon total energy deposited in HyCal
	double photon_energy = 0;
	photon_energy = pHandler->GetEnergy();

	for(iter=fecEventList->begin(); iter!=fecEventList->end(); ++iter)
	  {
            if( (*iter)->tag == 57633 ) // ti bank
              {
                vector<uint32_t> *ti_vec = (*iter)->getVector<uint32_t>();
                int ti_size = ti_vec->size();
		if(config.UseHyCalTimingCut == 1)
		  {
		    HyCalTimingCut = 0;
		    //HyCal timing cut
		    //G22  time slot 97  ; G1:112; G2:100; G6:113
		    //G21  time slot 115 ; G10:67; W12:92
		    //W31  time slot 116 ; W7:82 ; W8: 76
		    /***************************************************
		     * HyCal TDC Group 
		     * G1  112  :   G2  100 :  G3  109 :  G4  64
		     * G5  69   :   G6  113 :  G25 65  :  W6  66
		     * G10 67   :   W11  68 :  W18 70  :  W15 72
		     * W9  73   :   W20  74 :  W33 75  :  W8  76
		     * W27 77   :   W21  78 :   W2 80  :  W26 81
		     * W7  82   :   W1   83 :   W3 84  :  W14 85
		     * W24 88   :   G24  89 :   G20 90 :  G15 91
		     * W12 92   :   W36  93 :   W30 94 :  G11 96
		     * G22 97   :   W25  98 :   W13 99 :  G2  100
		     * W32 101  :   W28  104:   W4  105:  W10 106
		     * W16 107  :   W22  108:   G3  109:  G23 110
		     * G16 114  :  G21   115:   W31 116:  W19 117
		     * W23 120  :   W34  121:   W17 122:  W35 123
		     * W29 124  :    W5  125: Scin S1 126 : Scin S2 127
		     * ************************************************/
		    for(int i=0;i<ti_size;i++)
		      {
			if( (ti_vec->at(i) & 0xf8000000) !=0) continue;
			int tdc_ch = ( (ti_vec->at(i))>>19 )&0x7f;
			//if( (tdc_ch == 121)||(tdc_ch == 75)||(tdc_ch == 77) ||( tdc_ch == 104))
			//if( (tdc_ch == 67)||(tdc_ch == 92) )
			if( (tdc_ch == 123))
			  {
			    double tdc_value = (ti_vec->at(i)) & 0x7ffff; 
			    if(ntrigger<100) cout<<tdc_value<<endl;
			    if( (tdc_value>config.Hycal_Timing_Cut_Start) && (tdc_value<config.Hycal_Timing_Cut_End) ) 
			      {
				HyCalTimingCut = 1;  
				timing_test = tdc_value;
				break;
			      }
			  }
		      }
		  }

		// Scintillator Timing Cut ...
		// time slot 126
		// time slot 127
                if(config.UseScinTDCCut == 1)
		  {
		    convert = 0;
		    int TF126 = 0;
		    int TF127 = 0;
		    for(int i=0;i<ti_size;i++)
		      {
			if( (ti_vec->at(i) & 0xf8000000) !=0) continue;
			int tdc_ch = ( (ti_vec->at(i))>>19 )&0x7f;
		   
			//printf("0x%x \n", ti_vec->at(i));
			if(config.TDC_Channel == "126")
			  {
			    //cout<<" 126 cut..."<<endl;
			    if( (tdc_ch == 126)) 
			      { 
				double tdc_value = (ti_vec->at(i)) & 0x7ffff;
				//if( (tdc_value>7600) && (tdc_value<7800) )
				//if( (tdc_value>7000) && (tdc_value<10000) )
				//cout<<config.TDC_Start<<"  "<<config.TDC_End<<endl;
				if ( (tdc_value>config.TDC_Start) && (tdc_value<config.TDC_End) )
				  {
				    //cout<<"GEMPhys::ProcessAllEvents: tdc_value:  "<<tdc_value<<endl;
				    //hhTimeCorrelation->Fill(timing_test, tdc_value);
				    //hTimeDiff->Fill(tdc_value - timing_test);
				    convert = 1; 
				    break;
				  }
			      }
			  }

			if(config.TDC_Channel == "127")
			  {
			    //cout<<" 127 cut ..."<<endl;
			    if( (tdc_ch == 127)) 
			      { 
				double tdc_value = (ti_vec->at(i)) & 0x7ffff;
				//if( (tdc_value>7600) && (tdc_value<7800) )
				//if( (tdc_value>7000) && (tdc_value<10000) )
				if ( (tdc_value>config.TDC_Start) && (tdc_value<config.TDC_End) )
				  {
				    //cout<<"GEMPhys::ProcessAllEvents: tdc_value:  "<<tdc_value<<endl;
				    convert = 1; 
				    break;
				  }
			      }
			  }

			if(config.TDC_Channel == "126and127")
			  {
			    //cout<< " 126 and 127 cut ..."<<endl;
			    if( (tdc_ch == 126) ) TF126 = 1;
			    if( (tdc_ch == 127) ) TF127 = 1;
			    if( (TF126 == 1) && (TF127 == 1) ) 
			      { 
				double tdc_value = (ti_vec->at(i)) & 0x7ffff;
				//if( (tdc_value>7600) && (tdc_value<7800) )
				//if( (tdc_value>7000) && (tdc_value<10000) )
				if ( (tdc_value>config.TDC_Start) && (tdc_value<config.TDC_End) )
				  {
				    // cout<<"GEMPhys::ProcessAllEvents: tdc_value:  "<<tdc_value<<endl;
				    //hhTimeCorrelation->Fill(timing_test, tdc_value);
				    //hTimeDiff->Fill(tdc_value - timing_test);

				    convert = 1; 
				    break;
				  }
			      }
			  }

			if(config.TDC_Channel == "126or127")
			  {
			    //cout<< " 126 or 127 cut ..."<<endl;
			    if( (tdc_ch == 126) || (tdc_ch == 127) ) 
			      { 
				double tdc_value = (ti_vec->at(i)) & 0x7ffff;
				//if( (tdc_value>7600) && (tdc_value<7800) )
				//if( (tdc_value>7000) && (tdc_value<10000) )
				if ( (tdc_value>config.TDC_Start) && (tdc_value<config.TDC_End) )
				  {
				    //cout<<"GEMPhys::ProcessAllEvents: tdc_value:  "<<tdc_value<<endl;
				    hhTimeCorrelation->Fill(timing_test, tdc_value);
				    hTimeDiff->Fill(tdc_value - timing_test);

				    convert = 1; 
				    break;
				  }
			      }
			  }
		      }
		  }
              }
	  }

        if(ntrigger<10)
	  cout<<"HyCalTimingCut : "<<HyCalTimingCut
	      <<"  convert: "<<convert
	      <<"  photon energy: "<<photon_energy
	      <<endl;

	if( HyCalTimingCut == 1)
	  nHyCalEvents ++ ;
	if( convert == 1)
	  nScinEvents ++ ;

	if( (HyCalTimingCut == 1) && (convert == 1) && (photon_energy >= config.Hycal_Energy) ) 
	  { 
	    //cout<<"HyCal Energy: "<<config.Hycal_Energy<<endl;
	    if(ntrigger<10) cout<<"GEMPhys::ProcessEvents: Energy HyCal: "<<photon_energy<<endl;
	    nElectron+=1;
	  }
	else continue;
#endif

	for(iter=fecEventList->begin(); iter!=fecEventList->end(); ++iter)
	  {
#ifndef ZEROSUPPRESSION
	    if( ( (*iter)->tag == 57631) && 
	        ( ((int)((*iter)->num)==5) || ((int)((*iter)->num)==6) || 
		  ((int)((*iter)->num)==7) || ((int)((*iter)->num)==8) || 
		  ((int)((*iter)->num)==9) || ((int)((*iter)->num)==10)|| 
		  ((int)((*iter)->num)==11)|| ((int)((*iter)->num)==12) ) 
		)
	      {
		vector<uint32_t> *vec = (*iter)->getVector<uint32_t>();
		if(vec!=NULL  && vec->size()>0)
		  {
		    vSRSSingleEventData.reserve(vSRSSingleEventData.size() + vec->size() );
		    vSRSSingleEventData.insert(vSRSSingleEventData.end(), vec->begin(), vec->end() );
		  }
		else
		  {
		    cout<<"found NULL contents in fec.."<<endl;
		  }
	      }
#endif
#ifdef ZEROSUPPRESSION
            // zero suppression bank data
            if( ( (*iter)->tag == 57631) && ( ((int)(*iter)->num) == 99) )
	      {
		vector<uint32_t> *vec = (*iter)->getVector<uint32_t>();
		if(vec!=NULL  && vec->size()>0)
		  {
		    vSRSZeroEventData.reserve(vSRSZeroEventData.size() + vec->size() );
		    vSRSZeroEventData.insert(vSRSZeroEventData.end(), vec->begin(), vec->end() );
		  }
		else
		  {
		    cout<<"found NULL contents in fec.."<<endl;
		  }
	      }
#endif

	  }

	{	
	  ++entry;
	  if( (evtID != -1) && (entry > evtID) ) break; // process evtID many events
	  if(ntrigger<10) cout<<"GEMPhysHandelr::entry: "<<entry<< " nElectron: " << nElectron << " nScinEvents: "<<nScinEvents<<endl;
	  if(ntrigger<10) cout<<"SRS Event Size [uint_32]:"<<vSRSSingleEventData.size()<<endl;

#ifndef ZEROSUPPRESSION
	  if (vSRSSingleEventData.size() == 0 ) continue; // if no srs event found, go to next event
#endif

#ifdef ZEROSUPPRESSION
          if(vSRSZeroEventData.size() == 0 ) continue;
#endif
	  int size = vSRSSingleEventData.size();
	  uint32_t buf[size];
	  for(int i=0;i<size;i++)
	    {
	      buf[i]=vSRSSingleEventData[i];
	    }

          int zero_size = vSRSZeroEventData.size();
	  uint32_t buf_zero[zero_size];
	  for(int j=0;j<zero_size;j++)
	    {
	      buf_zero[j] = vSRSZeroEventData[j];
	    }

	  // do the work
#ifndef ZEROSUPPRESSION
	  GEMOnlineHitDecoder online_hit(buf, size, ped);
#endif
#ifdef ZEROSUPPRESSION
          GEMZeroHitDecoder online_hit(buf_zero, zero_size, ped);
#endif
	  // Fill histos
	  vector<GEMClusterStruct> gem1, gem2;
	  vector<GEMClusterStruct> hgem1, hgem2;
	  online_hit.GetClusterGEM(gem1, gem2);
	  online_hit.GetClusterHyCalCutMode(hgem1, hgem2);
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

	  //nb of clusters per event
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

	  // Calculate Detector Efficiency
	  if( (gem1.size()>0) || (gem2.size()>0) )
	    {
	      neff +=1.0;
	    }
	}
      }

    chan.close();

    /*
      TFile *f = new TFile("aaaa.root", "RECREATE");
      pHandler->GetTDCGroup("S1")->GetHist()->Write();
      pHandler->GetTDCGroup("S2")->GetHist()->Write();
      pHandler->GetTDCGroup("G21")->GetHist()->Write();
      pHandler->GetTDCGroup("G22")->GetHist()->Write();
      f->Close();
    */

  } catch (evioException e) {
    cerr <<endl <<e.toString() <<endl <<endl;
    exit(EXIT_FAILURE);
  }

  cout << "used time: " << timer.GetElapsedTime() << " ms" << endl;

  // claculate efficiency
  double eff = neff/ntrigger;
  cout<<"------------------------------------"<<endl;
  cout<<"Overall Detector Efficiency: "<<eff<<endl;
  cout<<"------------------------------------"<<endl;

#ifdef TDC_CHECK_EFFICIENCY
  double eff_real = (double)neff/nElectron;
  cout<<"GEM Hit: "<<neff<<endl;
  cout<<"TDC Events:  "<<nElectron<<endl;
  cout<<"Number of Events in Single File: "<<ntrigger<<endl;
  cout<<"------------------------------------"<<endl;
  cout<<"Overall Detector Efficiency: "<<eff_real<<endl;
  cout<<"------------------------------------"<<endl;
#endif

  return entry;
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

  vector<GEMClusterStruct> gem1, gem2;
  hit_decoder->GetClusterHyCalPlusMode(gem1, gem2);
  
  int s = HyCalGEMPosMatch(gem1, gem2, pHyCalHit);
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
  
  vector<GEMClusterStruct> gem1, gem2;
  hit_decoder->GetClusterHyCalPlusMode(gem1, gem2);

  GeometryMollerRing(gem1, gem2);
 
  int s = HyCalGEMPosMatch(gem1, gem2, pHyCalHit);

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
  
  vector<GEMClusterStruct> gem1, gem2;

  // GetCuster2DPositionBeamLine() returns a plane on GEM2, 
  // so the z coordinate is z_gem = z_gem2 = 5260mm
  hit_decoder->GetClusterHyCalPlusMode(gem1, gem2);

  //GeometryMollerRing(gem1, y1, x2, y2);
 
  int s = HyCalGEMPosMatch(gem1, gem2, pHyCalHit);

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

  vector<GEMClusterStruct> gem1, gem2;
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

int GEMPhysHandler::GEMHyCalPosMatch(int ngem, vector<GEMClusterStruct> &gem, vector<HyCalHit> *pHHit)
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
  
  vector<GEMClusterStruct> res_gem;

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

int GEMPhysHandler::HyCalGEMPosMatch( vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2, vector<HyCalHit> *pHHit)
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

  vector<GEMClusterStruct> res_gem1;
  vector<GEMClusterStruct> res_gem2;

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

void GEMPhysHandler::GeometryMollerRing(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct>& gem2)
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


