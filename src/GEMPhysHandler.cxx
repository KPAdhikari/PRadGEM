#include "GEMPhysHandler.h"
#include "GEMOnlineHitDecoder.h"
#include "GEMZeroHitDecoder.h"
#include "PRadDataHandler.h"
#include "PRadReconstructor.h"
#include "PRadTDCGroup.h"
#include "PRadBenchMark.h"
#include "PRadGEMTree.h"
#include <fstream>

using namespace std;
using namespace evio;

// GEM Efficiency Macros
#define TDC_CHECK_EFFICIENCY
#undef  TDC_CHECK_EFFICIENCY

#define ZEROSUPPRESSION
#undef ZEROSUPPRESSION

#define KEEP_MAX_TWO_CLUSTERS_AFTER_MATCH
#undef KEEP_MAX_TWO_CLUSTERS_AFTER_MATCH

#define PI 3.1415926

GEMPhysHandler::GEMPhysHandler()
{
    config = new GEMConfigure();
    config -> LoadConfigure();

    // run log
    string txt = config->phys_results_path;
    txt = txt + string(".run.log");
    outfile.open(txt.c_str(), std::ios::out);
    outfile.close();
    outfile.open(txt.c_str(), std::ios::in|std::ios::out);
    if( !outfile.is_open() )
	cout<<"Error: cannot open run_log file..."
	    <<endl;
    
    rst_tree = new PRadGEMTree();

    // tdc group log
    outfile<<"TDC group quantity: "
	<<config->TDC_Quan
	<<endl;
    for(int i=0;i<config->TDC_Quan;i++)
    {
	outfile<<"TDC Group: "
	    <<config->TDC[i]
	    <<endl;
    }

    InitTDCGroup();

    BookHistos();
    InitSetup();
    InitPedestal();
    InitVariables();
    InitHandler();
    InitIntermediateMollerCenterVariable();
}

void GEMPhysHandler::InitPedestal()
{
    string pedestal_file = config->GetLoadPedPath();
    ped = new GEMPedestal(pedestal_file);
    ped -> LoadPedestal( );
}

void GEMPhysHandler::InitSetup()
{
    // origin transfer
    O_Transfer = 253.2;
    OverlapStart = 231.2;
    Z_gem1 = 5300.+20.;//mm + chamber + screw cap thickness
    // need to consider:
    // ionization location depth, uncentiaty in z position

    Z_gem2 = 5260.+20.;//mm
    Z_hycal= 5600.;//mm
    Delta = 60.0;
    //Z_hycal = 5820;//mm
    beamEnergy = 2.147;//GeV
    //beamEnergy = 1.1;//GeV
    // 3% energy resolution 5 sigma
    beamEnergyCut = beamEnergy * 0.85 * 1000; //MeV
}

void GEMPhysHandler::InitVariables()
{
    vSRSSingleEventData.clear();
    vSRSZeroEventData.clear();
    FECs.clear();

    PRDMapping *mapping = PRDMapping::GetInstance();
    FECs = mapping->GetBankIDSet();
    //for debug
    outfile<<"FECs Found:  ";
    set<int>::iterator it;
    for(it=FECs.begin(); it!=FECs.end(); ++it)
    {
	outfile<<(*it)<<"  ";
    }
    outfile<<endl;

    totalEnergyDeposit = config->Hycal_Energy; //MeV

    nScinEvents = 0;
    nHyCalEvents = 0;

    // nElectron = nScinEvents && nHyCalEvents 
    nElectron = 0;

    // to be implemented...
    nElectron_126=0;
    nElectron_127=0;

    // to get how many events in total, 
    // this to compute photon conversion probability
    nTotalEvents = 0;  

    // # of events in gem, means gem events
    neff = 0;
    // # of gem events
    // after filter by hycal
    neff_after_match = 0;
    
    GEMMollerElectronQuantity = 0.;
    HyCalMollerElectronQuantity = 0.;
    GEMEpElectronQuantity = 0.;
    HyCalEpElectronQuantity = 0.;
    //sectorize gem
    for(int i=0;i<70;i++)
    {
        gem_moller_quantity[i] = 1.;
	hycal_moller_quantity[i] = 1.;
	gem_ep_quantity[i] = 0.;
	hycal_ep_quantity[i] = 0.;
    }
    //circular sectorize gem
    for(int i=0;i<50;i++)
    {
        r_sector[i] = 50. + 10.*i;
        gem_r_sec_ep_quantity[i] = 0.;
        hycal_r_sec_ep_quantity[i] = 0.;
    }
}

void GEMPhysHandler::InitHandler()
{
    pHandler = new PRadDataHandler();
    pHandler->ReadTDCList("/home/xbai/w/pRad/source/PRadDecoder/config/tdc_group_list.txt");
    pHandler->ReadChannelList("/home/xbai/w/pRad/source/PRadDecoder/config/module_list.txt");
    pHandler->BuildChannelMap();
    pHandler->ReadPedestalFile("/home/xbai/w/pRad/source/PRadDecoder/config/pedestal.dat");
    pHandler->ReadCalibrationFile("/home/xbai/w/pRad/source/PRadDecoder/config/calibration.txt");
    //pHandler->EnableReconstruction();

    reconstruct = new PRadReconstructor();
    reconstruct->SetHandler(pHandler);
}

void GEMPhysHandler::InitIntermediateMollerCenterVariable()
{
    //compute intersection points
    px1 = 0x270F;py1=0x270F;px2=0x270F;py2=0x270F;
    cx1 = 0x270F;cy1=0x270F;cx2=0x270F;cy2=0x270F;

    px1_c = 0x270F;py1_c=0x270F;px2_c=0x270F;py2_c=0x270F;
    cx1_c = 0x270F;cy1_c=0x270F;cx2_c=0x270F;cy2_c=0x270F;
}

GEMPhysHandler::~GEMPhysHandler()
{
    outfile.close();
    ped->Delete();
}

void GEMPhysHandler::ProcessAllFiles()
{
    int nFile = config->nFile;
    outfile<<"GEMPhysHandler::ProcessAllFiles: # of Files to analyze: "
	<<nFile
	<<endl;

    for(int i=0;i<nFile;i++)
    {
	outfile<<config->fileList[i]<<endl;
    }

    outfile<<"TDC Cut: "<<config->TDC_Channel
	<<", Start TIME:  "<<config->TDC_Start
	<<", END TIME:  "<<config->TDC_End
	<<", HyCal Energy Cut:  "<<config->Hycal_Energy
	<<endl;

    // to compute photon conversion rate
    nTotalEvents = 0;
    neff = 0;
    neff_after_match = 0;
    nElectron = 0;

    if(config->fileList[0].find("evio.0")!= string::npos)
    {
	pHandler->InitializeByData(config->fileList[0].c_str());
    }

    cout<<"cache file name list..."<<endl;
    string file_list[nFile];
    for(int i=0;i!=nFile;i++)
    {
        file_list[i] = config->fileList[i];
    }

    for(int i=1;i!=nFile;++i)
    {
        string nam = file_list[i];
	filename = nam;
	ProcessAllEvents(-1);
	//ProcessAllEvents(4000);
    }

    // Save Histos
    SavePhysResults();

    outfile<<"<><><><><><><>"<<endl;
    outfile<<"Overall Detector Efficiency in production runs:"
           <<endl;
    outfile<<"Moller Events Efficiency: "
           <<GEMMollerElectronQuantity/HyCalMollerElectronQuantity
	   <<endl;
    outfile<<"e-p Events Efficiency: "
           <<GEMEpElectronQuantity<<" / "
	   <<HyCalEpElectronQuantity<<" = "
           <<GEMEpElectronQuantity/HyCalEpElectronQuantity
	   <<endl;
    WriteSectorEff();

#ifdef TDC_CHECK_EFFICIENCY
    outfile<<endl<<endl;
    double eff_real = (double)neff/nElectron;
    double eff_real_after_match = (double)neff_after_match/nElectron;
    outfile<<"GEM Hit: "<<neff<<endl;
    outfile<<"GEM Hit after match: "<<neff_after_match<<endl;
    outfile<<"TDC:  "<<nElectron<<endl;
    outfile<<"Total Number of Events in "<<nFile
	<<" Files: "<<nTotalEvents
	<<endl;
    outfile<<"------------------------------------"<<endl;
    outfile<<"Overall Detector Efficiency: "
           <<eff_real
	   <<endl;
    outfile<<"Overall Detector Efficiency after match: "
           <<eff_real_after_match
	   <<endl;
    outfile<<"------------------------------------"<<endl;
#endif
}

int GEMPhysHandler::ProcessAllEvents(int evtID )
{
    outfile<<"Process File:  "<<filename<<endl;
    cout<<"Process File:  "<<filename<<endl;

    int entry = 0;
    double neff_previous = neff;
    double neff_previous_after_match = neff_after_match;
    double nElectron_previous = nElectron;
    double ntrigger_current_file = 0.0; // for detector efficiency

    //for detector efficiency in producton run
    double hycal_moller_quantity_previous = HyCalMollerElectronQuantity;
    double hycal_ep_quantity_previous = HyCalEpElectronQuantity;
    double gem_moller_quantity_previous = GEMMollerElectronQuantity;
    double gem_ep_quantity_previous = GEMEpElectronQuantity;

    PRadBenchMark timer;

    try{
	evioFileChannel chan(filename.c_str(), "r");
	chan.open();
	while(chan.read())
	{
	    ntrigger_current_file +=1.0;
	    nTotalEvents += 1.0;

	    pHandler->Decode(chan.getBuffer());
	    pHyCalHit =& reconstruct->CoarseHyCalReconstruct(pHandler->GetEventCount() - 1); 

	    hQuantityOfClustersHyCal->Fill(pHyCalHit->size());
	    for(int i=0;i<pHyCalHit->size(); i++)
	    {
		hhHyCalClusterMap->Fill(pHyCalHit->at(i).x, pHyCalHit->at(i).y);

		hHyCalEnergy->Fill(pHyCalHit->at(i).E);
		if(pHyCalHit->size() == 1) 
		    hHyCalEnergyEp->Fill(pHyCalHit->at(i).E);
		if(pHyCalHit->size() == 2) 
		{
		    hHyCalEnergyMoller->Fill(pHyCalHit->at(i).E);
		}
	    }

	    vSRSSingleEventData.clear();
	    vSRSZeroEventData.clear();

	    evioDOMTree event(chan);
	    evioDOMNodeListP fecEventList = event.getNodeList( isLeaf() );
	    //outfile<<"total number of all banks: "<<fecEventList->size()<<endl;

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
		    GetCutFlags(ntrigger_current_file, ti_vec, convert, HyCalTimingCut);
		}
	    }

	    if(ntrigger_current_file<10)
		outfile<<"HyCalTimingCut : "<<HyCalTimingCut
		    <<"  convert: "<<convert
		    <<"  photon energy: "<<photon_energy
		    <<endl;

	    if( HyCalTimingCut == 1)
		nHyCalEvents ++ ;
	    if( convert == 1)
		nScinEvents ++ ;

	    if( (HyCalTimingCut == 1) && (convert == 1) && (photon_energy >= config->Hycal_Energy) ) 
	    { 
		//outfile<<"HyCal Energy: "<<config->Hycal_Energy<<endl;
		if(ntrigger_current_file<10) 
		    outfile<<"GEMPhys::ProcessEvents: Energy HyCal: "
			<<photon_energy
			<<endl;
		nElectron+=1;
	    }
	    else 
		continue;
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

	    ++entry;
	    if( (evtID != -1) && (entry > evtID) ) 
		break; // process evtID many events

	    if(ntrigger_current_file<10) 
		outfile<<"GEMPhysHandelr::entry: "<<entry
		    << " nElectron: " << nElectron 
		    << " nScinEvents: "<<nScinEvents
		    <<endl;
	    if(ntrigger_current_file<10) 
		outfile<<"SRS Event Size [uint_32]:"<<vSRSSingleEventData.size()
		    <<endl;

#ifndef ZEROSUPPRESSION
	    if (vSRSSingleEventData.size() == 0 ) 
		continue; 
#endif

#ifdef ZEROSUPPRESSION
	    if(vSRSZeroEventData.size() == 0 ) 
		continue;
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

#ifndef ZEROSUPPRESSION
	    GEMOnlineHitDecoder online_hit(buf, size, ped);
#endif
#ifdef ZEROSUPPRESSION
	    GEMZeroHitDecoder online_hit(buf_zero, zero_size, ped);
#endif

	    // Fill histos
	    ComputeGEMOffsets(&online_hit);
            EvalMatchMech(&online_hit);
	    ProcessEp(&online_hit);
	    ProcessMoller(&online_hit);
	    ProcessMollerAfterCorrection(&online_hit);
	    GetGEMClusterMapHyCalCoor(&online_hit);
	    CharactorizeGEM(&online_hit);
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

    outfile << "used time: " << timer.GetElapsedTime() << " ms" << endl;
    cout << "used time: " << timer.GetElapsedTime() << " ms" << endl;

    // claculate efficiency
    double eff = (double)(neff-neff_previous)/ntrigger_current_file;
    double eff_after_match = (double)(neff_after_match-neff_previous_after_match)/ntrigger_current_file;
    outfile<<"photon conversion rate in single file: "<<eff<<endl;
    outfile<<"photon conversion rate in single file after match: "<<eff_after_match<<endl;
    outfile<<"<><><><><><><><>"
           <<endl;
    outfile<<"Single File GEM Efficiency in producton runs: "
           <<endl;
    outfile<<"moller events efficiency: "
           <<(GEMMollerElectronQuantity - gem_moller_quantity_previous) / (HyCalMollerElectronQuantity - hycal_moller_quantity_previous)
	   <<endl;
    outfile<<"ep events efficiency: "
           <<(GEMEpElectronQuantity - gem_ep_quantity_previous) / (HyCalEpElectronQuantity - hycal_ep_quantity_previous)
	   <<endl;


#ifdef TDC_CHECK_EFFICIENCY
    double eff_real = (double)(neff-neff_previous)/(nElectron-nElectron_previous);
    double eff_real_after_match=(double)(neff_after_match-neff_previous_after_match)/(nElectron-nElectron_previous);
    outfile<<"GEM Hit: "
           <<neff-neff_previous
	   <<endl;
    outfile<<"GEM Hit after match: "
           <<neff_after_match-neff_previous_after_match
	   <<endl;
    outfile<<"Scin & HyCal TDC cut Events:  "
           <<nElectron-nElectron_previous
	   <<endl;
    outfile<<"Number of Events in Single File: "
           <<ntrigger_current_file
	   <<endl;
    outfile<<"Detector Efficiency in single file : "
           <<eff_real
	   <<endl<<endl;
    outfile<<"Detector Efficiency in single file after match : "
           <<eff_real_after_match
	   <<endl<<endl;
#endif
    return entry;
}

int GEMPhysHandler::GetTDCGroup( const string &st)
{
    return TDC_Map[st]; 
}

void GEMPhysHandler::InitTDCGroup()
{
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

    TDC_Map[string("G1")] = 112;TDC_Map[string("G2")] = 100;TDC_Map[string("G3")] = 109;TDC_Map[string("G4")] = 64;
    TDC_Map[string("G5")] = 69;TDC_Map[string("G6")] = 113;TDC_Map[string("G25")] = 65;TDC_Map[string("W6")] = 66;
    TDC_Map[string("G10")] = 67;TDC_Map[string("W11")] = 68;TDC_Map[string("W18")] = 70;TDC_Map[string("W15")] = 72;
    TDC_Map[string("W9")] = 73;TDC_Map[string("W20")] = 74;TDC_Map[string("W33")] = 75;TDC_Map[string("W8")] = 76;
    TDC_Map[string("W27")] = 77;TDC_Map[string("W21")] = 78;TDC_Map[string("W2")] = 80;TDC_Map[string("W26")] = 81;
    TDC_Map[string("W7")] = 82;TDC_Map[string("W1")] = 83;TDC_Map[string("W3")] = 84;TDC_Map[string("W14")] = 85;
    TDC_Map[string("W24")] = 88;TDC_Map[string("G24")] = 89;TDC_Map[string("G20")] = 90;TDC_Map[string("G15")] = 91;
    TDC_Map[string("W12")] = 92;TDC_Map[string("W36")] = 93;TDC_Map[string("W30")] = 94;TDC_Map[string("G11")] = 96;
    TDC_Map[string("G22")] = 97;TDC_Map[string("W25")] = 98;TDC_Map[string("W13")] = 99;TDC_Map[string("G2")] = 100;
    TDC_Map[string("W32")] = 101;TDC_Map[string("W28")] = 104;TDC_Map[string("W4")] = 105;TDC_Map[string("W10")] = 106;
    TDC_Map[string("W16")] = 107;TDC_Map[string("W22")] = 108;TDC_Map[string("G3")] = 109;TDC_Map[string("G23")] = 110;
    TDC_Map[string("G16")] = 114;TDC_Map[string("G21")] = 115;TDC_Map[string("W31")] = 116;TDC_Map[string("W19")] = 117;
    TDC_Map[string("W23")] = 120;TDC_Map[string("W34")] = 121;TDC_Map[string("W17")] = 122;TDC_Map[string("W35")] = 123;
    TDC_Map[string("W29")] = 124;TDC_Map[string("W5")] = 125;TDC_Map[string("S1")] = 126;TDC_Map[string("S2")] = 127;
}  

void GEMPhysHandler::GetCutFlags(int nth_event, vector<uint32_t> *ti_vec, int & convert, int & HyCalTimingCut)
{
    int ti_size = ti_vec->size();

    if(config->UseHyCalTimingCut == 1)
    {

	HyCalTimingCut = 0;
	for(int i=0;i<ti_size;i++)
	{
	    if( (ti_vec->at(i) & 0xf8000000) !=0) continue;
	    int tdc_ch = ( (ti_vec->at(i))>>19 )&0x7f;

	    int right_tdc_channel = 0;
	    for(int i=0;i<config->TDC_Quan;i++)
	    {
		if(tdc_ch == GetTDCGroup(config->TDC[i]) )
		    right_tdc_channel = 1;
	    }

	    if( right_tdc_channel == 1)
	    {
		double tdc_value = (ti_vec->at(i)) & 0x7ffff; 

		if(nth_event<100) 
		    outfile<<tdc_value<<endl;

		if( (tdc_value>config->Hycal_Timing_Cut_Start) && (tdc_value<config->Hycal_Timing_Cut_End) ) 
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
    if(config->UseScinTDCCut == 1)
    {
	convert = 0;
	int TF126 = 0;
	int TF127 = 0;
	for(int i=0;i<ti_size;i++)
	{
	    if( (ti_vec->at(i) & 0xf8000000) !=0) continue;
	    int tdc_ch = ( (ti_vec->at(i))>>19 )&0x7f;

	    //printf("0x%x \n", ti_vec->at(i));
	    if(config->TDC_Channel == "126")
	    {
		//outfile<<" 126 cut..."<<endl;
		if( (tdc_ch == 126)) 
		{ 
		    double tdc_value = (ti_vec->at(i)) & 0x7ffff;
		    //if( (tdc_value>7600) && (tdc_value<7800) )
		    //if( (tdc_value>7000) && (tdc_value<10000) )
		    //outfile<<config->TDC_Start<<"  "<<config->TDC_End<<endl;
		    if ( (tdc_value>config->TDC_Start) && (tdc_value<config->TDC_End) )
		    {
			//outfile<<"GEMPhys::ProcessAllEvents: tdc_value:  "<<tdc_value<<endl;
			//hhTimeCorrelation->Fill(timing_test, tdc_value);
			//hTimeDiff->Fill(tdc_value - timing_test);
			convert = 1; 
			break;
		    }
		}
	    }

	    if(config->TDC_Channel == "127")
	    {
		//outfile<<" 127 cut ..."<<endl;
		if( (tdc_ch == 127)) 
		{ 
		    double tdc_value = (ti_vec->at(i)) & 0x7ffff;
		    //if( (tdc_value>7600) && (tdc_value<7800) )
		    //if( (tdc_value>7000) && (tdc_value<10000) )
		    if ( (tdc_value>config->TDC_Start) && (tdc_value<config->TDC_End) )
		    {
			//outfile<<"GEMPhys::ProcessAllEvents: tdc_value:  "<<tdc_value<<endl;
			convert = 1; 
			break;
		    }
		}
	    }

	    if(config->TDC_Channel == "126and127")
	    {
		//outfile<< " 126 and 127 cut ..."<<endl;
		if( (tdc_ch == 126) ) TF126 = 1;
		if( (tdc_ch == 127) ) TF127 = 1;
		if( (TF126 == 1) && (TF127 == 1) ) 
		{ 
		    double tdc_value = (ti_vec->at(i)) & 0x7ffff;
		    //if( (tdc_value>7600) && (tdc_value<7800) )
		    //if( (tdc_value>7000) && (tdc_value<10000) )
		    if ( (tdc_value>config->TDC_Start) && (tdc_value<config->TDC_End) )
		    {
			// outfile<<"GEMPhys::ProcessAllEvents: tdc_value:  "<<tdc_value<<endl;
			//hhTimeCorrelation->Fill(timing_test, tdc_value);
			//hTimeDiff->Fill(tdc_value - timing_test);

			convert = 1; 
			break;
		    }
		}
	    }

	    if(config->TDC_Channel == "126or127")
	    {
		//outfile<< " 126 or 127 cut ..."<<endl;
		if( (tdc_ch == 126) || (tdc_ch == 127) ) 
		{ 
		    double tdc_value = (ti_vec->at(i)) & 0x7ffff;
		    //if( (tdc_value>7600) && (tdc_value<7800) )
		    //if( (tdc_value>7000) && (tdc_value<10000) )
		    if ( (tdc_value>config->TDC_Start) && (tdc_value<config->TDC_End) )
		    {
			//outfile<<"GEMPhys::ProcessAllEvents: tdc_value:  "<<tdc_value<<endl;
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

template<class T> void GEMPhysHandler::ProcessEp(T * hit_decoder)
{
    // theta distribution
    double z_gem1 = Z_gem1;
    double z_gem2 = Z_gem2;
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
	if(gem1[0].energy < totalEnergyDeposit ) 
            return;
	gem1[0].x -= xoffset_beam;
	gem1[0].y -= yoffset_beam;

	theta = TMath::Sqrt((gem1[0].x)*(gem1[0].x) + gem1[0].y*gem1[0].y) / z_gem1;
	theta = TMath::ATan(theta);
	double theta_d = theta*180.0/PI;
	
        hThetaDistribution->Fill(theta_d);

	hThetaDistributionEp->Fill(theta_d);
	hhEnergyVsAngle->Fill(theta_d, gem1[0].energy);
	hhEnergyVsAngleEp->Fill(theta_d, gem1[0].energy);
	if( gem1[0].energy >= beamEnergyCut )
	{
	    hhEnergyVsAngleEnergyCut->Fill(theta_d, gem1[0].energy);
	    hhEnergyVsAngleEpEnergyCut->Fill(theta_d, gem1[0].energy);
	}

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
    else if( (gem2.size() == 1)&&(gem1.size() == 0))
    {
	//Ep energy cut
	if(gem2[0].energy < totalEnergyDeposit) return;

	gem2[0].x -= xoffset_beam;
	gem2[0].y -= yoffset_beam;

	theta = TMath::Sqrt( (gem2[0].x)*(gem2[0].x) + gem2[0].y*gem2[0].y) / z_gem2;
	theta = TMath::ATan(theta);
	double theta_d = theta*180.0/PI;

	hThetaDistribution->Fill(theta_d);

	hThetaDistributionEp->Fill(theta_d);
	hhEnergyVsAngle->Fill(theta_d, gem2[0].energy);
	hhEnergyVsAngleEp->Fill(theta_d, gem2[0].energy);
	if( gem2[0].energy >= beamEnergyCut )
	{
	    hhEnergyVsAngleEnergyCut->Fill(theta_d, gem2[0].energy);
	    hhEnergyVsAngleEpEnergyCut->Fill(theta_d, gem2[0].energy);
	}


	//q suqare
	top = 4.0*beamEnergy*beamEnergy*TMath::Sin(theta/2)*TMath::Sin(theta/2);
	bottom = 1+ (2*beamEnergy/0.938)*TMath::Sin(theta/2)*TMath::Sin(theta/2);
	q_square = top/bottom;
	hQSquareEp->Fill(q_square); 
	hhQSquareScattAngleEp -> Fill(theta_d, q_square);

	if((theta_d > 0.7) && (theta_d<0.8) ) 
	    hQSquareEp1->Fill(q_square);

	if((theta_d > 1.0)&& (theta_d<1.1) ) 
	    hQSquareEp2->Fill(q_square);

	if((theta_d > 1.5)&& (theta_d<1.6) ) 
	    hQSquareEp3->Fill(q_square);

	if((theta_d > 2.0)&& (theta_d<2.1) ) 
	    hQSquareEp4->Fill(q_square);

	if(theta_d > 2.2) 
	    hQSquareEp5->Fill(q_square);

    }
}

template<class T> void GEMPhysHandler::ProcessMoller(T * hit_decoder)
{
    // theta distribution
    double z_gem1 = Z_gem1;
    double z_gem2 = Z_gem2;
    double theta = 0;
    double q_square = 0;
    double top = 0.0;
    double bottom = 0.0; //q_square = top/bottom
    double e0 = beamEnergy; //GeV

    vector<GEMClusterStruct> gem1, gem2;
    hit_decoder->GetClusterHyCalPlusMode(gem1, gem2);

    GeometryMollerRing(gem1, gem2);

    int s = HyCalGEMPosMatch(gem1, gem2, pHyCalHit);

    if( s == 0 ) 
        return ;

    vector<GEMClusterStruct> gem;
    for(int i=0;i<gem1.size();i++)
    {
        gem1[i].SetZ( (float) z_gem1);
        gem.push_back(gem1[i]);
    }
    for(int i=0;i<gem2.size();i++)
    {
        gem2[i].SetZ( (float) z_gem2);
        gem.push_back(gem2[i]);
    }

    if(gem.size() != 2)
        return;

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
    if( (gem.size() == 2) ) 
    { 
	if(  ( (gem[0].x)*(gem[1].x) < 0 ) && ( gem[0].y*gem[1].y < 0 ) ) 
	{
	    if( (gem[0].energy+gem[1].energy) < totalEnergyDeposit) 
	        return;

	    double temp = 0;
	    double slope1 = 0;
	    double slope2 = 0;
	    //1st electron
	    theta = TMath::Sqrt((gem[0].x)*(gem[0].x) + gem[0].y*gem[0].y) / gem[0].z;
	    theta = TMath::ATan(theta);
	    hThetaDistribution->Fill(theta*180.0/PI);

	    hhEnergyVsAngle->Fill(theta*180.0/PI, gem[0].energy);
	    hhEnergyVsAngleMoller->Fill(theta*180.0/PI, gem[0].energy);

	    slope1 = TMath::Abs( gem[0].y/gem[0].x );
	    slope1 = TMath::ATan(slope1); 
	    assert( slope1 < PI/2);
	    if( (gem[0].x<0) && (gem[0].y>0) ) 
		slope1 = PI - slope1; 
	    else if( (gem[0].x<0) && (gem[0].y<0) ) 
		slope1 = slope1+PI; 
	    else if( (gem[0].x>0) && (gem[0].y <0)) 
		slope1 = 2*PI - slope1 ;

	    temp = theta;

	    //2nd electron
	    theta = TMath::Sqrt((gem[1].x)*(gem[1].x) + gem[1].y*gem[1].y) / gem[1].z;
	    theta = TMath::ATan(theta);
	    hThetaDistribution->Fill(theta*180.0/PI);

	    hhEnergyVsAngle->Fill(theta*180.0/PI, gem[1].energy);
	    hhEnergyVsAngleMoller->Fill(theta*180.0/PI, gem[1].energy);

	    if( gem[0].energy + gem[1].energy > beamEnergyCut)
	    {
		hhEnergyVsAngleEnergyCut->Fill(temp*180.0/PI, gem[0].energy);
	        hhEnergyVsAngleMollerEnergyCut->Fill(temp*180.0/PI, gem[0].energy);
        	hhEnergyVsAngleEnergyCut->Fill(theta*180.0/PI, gem[1].energy);
	        hhEnergyVsAngleMollerEnergyCut->Fill(theta*180.0/PI, gem[1].energy);
	    }

	    slope2 = TMath::Abs( gem[1].y / (gem[1].x) );
	    slope2 = TMath::ATan(slope2);
	    if( (gem[1].x<0) && (gem[1].y>0) ) 
		slope1 = PI - slope1; 
	    else if( (gem[1].x<0) && (gem[1].y<0) ) 
		slope1 = slope1+PI; 
	    else if( (gem[1].x>0) && (gem[1].y <0)) 
		slope1 = 2*PI - slope1;

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
	    if( ( (temp>ThetaSmall) && (temp<ThetaLarge)) || ( (theta>ThetaSmall)&&(theta<ThetaLarge)   ) )
	    {
		hhMoller2DRing->Fill(gem[0].x, gem[0].y);
		hhMoller2DRing->Fill(gem[1].x, gem[1].y);
	    }
	    if( ( (temp>ThetaSmall2) && (temp<ThetaLarge2)) || ( (theta>ThetaSmall2)&&(theta<ThetaLarge2)   ) )
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
		    cx2 = gem[1].x;
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
			double xi = (b1*c2 - b2*c1)/D; //outfile<<xi<<endl;
			double yi = (a2*c1 - a1*c2)/D; //outfile<<yi<<endl;
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
}

template<class T> void GEMPhysHandler::ProcessMollerAfterCorrection(T * hit_decoder)
{
    // theta distribution
    double z_gem1 = Z_gem1;
    double z_gem2 = Z_gem2;
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

    vector<GEMClusterStruct> gem;
    for(int i=0;i<gem1.size();i++)
    {
        gem1[i].SetZ(z_gem1);
	gem.push_back(gem1[i]);
    }
    for( int i=0;i<gem2.size();i++)
    {
        gem2[i].SetZ(z_gem2);
	gem.push_back(gem2[i]);
    }
     
    if(gem.size() != 2)
        return;

    // origin offset from beam line
    // z coordinate: GEM2 : 5260mm
    double xoffset_beam = 1.631;
    double yoffset_beam = 0.366;


    // Moller events selection
    // require them in different quadrants, very rough

    if( gem.size() == 2 ) 
    { 
	if(  ( (gem[0]).x*(gem[1]).x < 0 ) && ( gem[0].y*gem[1].y < 0 ) ) 
	{
	    //beamline correction
	    gem[0].x -= xoffset_beam;
	    gem[1].x -= xoffset_beam;
	    gem[0].y -= yoffset_beam;
	    gem[1].y -= yoffset_beam;

	    if( (gem[0].energy+gem[1].energy) < totalEnergyDeposit) 
	        return;

	    double temp = 0;
	    double slope1 = 0;
	    double slope2 = 0;

	    //1st electron
	    theta = TMath::Sqrt((gem[0].x)*(gem[0].x) + gem[0].y*gem[0].y) / gem[0].z;
	    theta = TMath::ATan(theta);
	    hhEnergyVsAngleBeamLineCorrection->Fill(theta*180.0/PI, gem[0].energy);
	    hhEnergyVsAngleMollerBeamLineCorrection->Fill(theta*180.0/PI, gem[0].energy);

	    slope1 = TMath::Abs( gem[0].y/gem[0].x );
	    slope1 = TMath::ATan(slope1); 
	    assert( slope1 < PI/2);
	    if( (gem[0].x<0) && (gem[0].y>0) ) 
	        slope1 = PI - slope1; 
	    else if( (gem[0].x<0) && (gem[0].y<0) ) 
	        slope1 = slope1+PI; 
	    else if( (gem[0].x>0) && (gem[0].y <0)) 
	        slope1 = 2*PI - slope1 ;
	    temp = theta;

	    //2nd electron
	    theta = TMath::Sqrt((gem[1].x)*(gem[1].x) + gem[1].y*gem[1].y) / gem[1].z;
	    theta = TMath::ATan(theta);
	    hhEnergyVsAngleBeamLineCorrection->Fill(theta*180.0/PI, gem[1].energy);
	    hhEnergyVsAngleMollerBeamLineCorrection->Fill(theta*180.0/PI, gem[1].energy);

	    slope2 = TMath::Abs( gem[1].y / (gem[1].x) );
	    slope2 = TMath::ATan(slope2);
	    if( (gem[1].x<0) && (gem[1].y>0) ) 
		slope1 = PI - slope1; 
	    else if( (gem[1].x<0) && (gem[1].y<0) ) 
		slope1 = slope1+PI; 
	    else if( (gem[1].x>0) && (gem[1].y<0)) 
		slope1 = 2*PI - slope1 ;
	    
	    // moller angular resolution
	    rst_tree -> Clear();
            if( (gem[0].energy>485.*beamEnergy) && (gem[0].energy<515.*beamEnergy) && 
	        (gem[1].energy>485.*beamEnergy) && (gem[1].energy<515.*beamEnergy) )
	    {
	        rst_tree -> angular_resolution_from_moller1 = temp*180./PI - MollerAngleFromEnergy(gem[0].energy);
	        rst_tree -> angular_resolution_from_moller2 = theta*180./PI - MollerAngleFromEnergy(gem[1].energy);
		rst_tree -> angle_moller1 = temp*180./PI;
		rst_tree -> angle_moller2 = temp*180./PI;
		rst_tree -> symm_dx = gem[0].x + gem[1].x;
		rst_tree -> symm_dy = gem[0].y + gem[1].y;
		rst_tree -> symm_coplanarity = (slope1-slope2)*180.0/PI - 180.0;
	    }

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
		    cx1_c = gem[0].x;
		    cx2_c = gem[1].x;
		    cy1_c = gem[0].y;
		    cy2_c = gem[1].y;
		    double a1 = py1_c - py2_c;
		    double b1 = px2_c - px1_c;
		    double c1 = px1_c * py2_c - px2_c * py1_c;

		    double a2 = cy1_c - cy2_c;
		    double b2 = cx2_c - cx1_c;
		    double c2 = cx1_c * cy2_c - cx2_c * cy1_c;

		    double D = a1*b2 - a2*b1;

		    if(D != 0)
		    {
			double xi = (b1*c2 - b2*c1)/D; //outfile<<xi<<endl;
			double yi = (a2*c1 - a1*c2)/D; //outfile<<yi<<endl;
			hXOffsetFromMollerAfterCorrection->Fill(xi);
			hYOffsetFromMollerAfterCorrection->Fill(yi);
			hhMollerCenterAfterCorrection->Fill(xi, yi);
		    }

		}
		else
		{
		    px1_c = cx1_c = gem[0].x;
		    py1_c = cy1_c = gem[0].y;
		    px2_c = cx2_c = gem[1].x;
		    py2_c = cy2_c = gem[1].y;
		}
	    }
	} 
    }
}

template<class T> void GEMPhysHandler::GetGEMClusterMapHyCalCoor(T * online_hit)
{
    vector<GEMClusterStruct> hgem1, hgem2;
    online_hit->GetClusterHyCalCutMode(hgem1, hgem2);

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


template<class T> void GEMPhysHandler::ComputeGEMOffsets(T * online_hit)
{
    vector<GEMClusterStruct> gem1, gem2;
    online_hit->GetClusterGEM(gem1, gem2);

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
		double z_gem1 = Z_gem1;
		double z_gem2 = Z_gem2;
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

    if( gem1.size()>0 || gem2.size()>0)
    {
	online_hit->FillHistos(hNbClusterPerPlaneX, hNbClusterPerPlaneY, hClusterDistX, hClusterDistY);

	// for detector efficiency
	neff += 1.0;
    }

    // efficiency after match
    gem1.clear();
    gem2.clear();
    online_hit->GetClusterHyCalCutMode(gem1, gem2);
    GEMHyCalPosMatch(0, gem1, pHyCalHit);
    GEMHyCalPosMatch(1, gem2, pHyCalHit);
    if( gem1.size() + gem2.size() > 0)
    {
        neff_after_match += 1.0;
	for(int i=0;i<gem1.size();i++)
	    hhGEMClusterMapAfterMatch->Fill(gem1[i].x, gem1[i].y);
	for(int i=0;i<gem2.size();i++)
	    hhGEMClusterMapAfterMatch->Fill(gem2[i].x, gem2[i].y);
    }

}

template<class T> void GEMPhysHandler::CharactorizeGEM(T * hit_decoder)
{
    // theta distribution
    double z_gem1 = Z_gem1;
    double z_gem2 = Z_gem2;
    double theta = 0;
    double q_square = 0;
    double top = 0.0;
    double bottom = 0.0; //q_square = top/bottom
    double e0 = beamEnergy; //GeV

    vector<GEMClusterStruct> gem1, gem2;
    hit_decoder->GetClusterGEM(gem1, gem2);
   
    rst_tree -> Fill(nTotalEvents, gem1, gem2, *pHyCalHit);

    if(gem1.size() > 0)
    {
	for(int i=0;i<gem1.size();i++)
	{
	    //theta = TMath::Sqrt( (gem1[i].x-O_Transfer)*(gem1[i].x-O_Transfer) + gem1[i].y*gem1[i].y) / z_gem1;
	    //theta = TMath::ATan(theta); 
	    //hThetaDistribution->Fill(theta*180.0/PI);

	    hhClusterDist[0] -> Fill( gem1[i].x, gem1[i].y);
	    hhClusterDist_pRad[0]->Fill( gem1[i].x, (-1.0)*gem1[i].y );
	    hhChargeRatio[0] -> Fill(gem1[i].x_charge, gem1[i].y_charge);
	    hChargeRatio[0] -> Fill( gem1[i].x_charge/gem1[i].y_charge);
	    hClusterADCDistXSide[0]->Fill( gem1[i].x_charge);
	    hClusterADCDistYSide[0]->Fill( gem1[i].y_charge);
	    if( gem1[i].x_charge<0)
		outfile<<"x1 negative charge???"<<endl;;
	    if( gem1[i].y_charge<0)
		outfile<<"y1 negative charge???"<<endl;;

	    // save tree 
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
	    if( gem2[i].x_charge<0)
		outfile<<"x2 negative charge???"<<endl;
	    if( gem2[i].y_charge<0)
		outfile<<"y2 negative charge???"<<endl;

             // save tree
	}
    }
}

void GEMPhysHandler::SavePhysResults()
{
    // save root 
    rst_tree -> Save();
    WriteHistos();

}

int GEMPhysHandler::GEMHyCalPosMatch(int ngem, vector<GEMClusterStruct> &gem, vector<HyCalHit> *pHHit)
{
    // filter GEM clusters based on HyCal clusters

    if(gem.size() == 0)
    {
	//pHHit->clear();
	return 0;
    }
    if(pHHit->size() == 0) 
    {
	gem.clear();
	return 0;
    }

    double z_gem1 = Z_gem1;
    double z_gem2 = Z_gem2;
    double z_gem = 0;
    double z_hycal = Z_hycal;
    double res = Delta; // a larger range, 60mm

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
	    double dr = TMath::Sqrt( (x_hycal-pHHit->at(j).x)*(x_hycal-pHHit->at(j).x) + 
	                             (y_hycal-pHHit->at(j).y)*(y_hycal-pHHit->at(j).y)  );

            hROffsetGEMHyCal->Fill(dr);
	    if( (dr - res)<=0.)
	    {
		gem[i].energy = pHHit->at(j).E;
		res_gem.push_back(gem[i]);
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

int GEMPhysHandler::HyCalGEMPosMatch( vector<GEMClusterStruct> &gem1, 
                                      vector<GEMClusterStruct> &gem2, 
				      vector<HyCalHit> *pHHit)
{
    if( (gem1.size() == 0) && (gem2.size() == 0))
    {
	return 0;
    }
    if(pHHit->size() == 0) 
    {
	gem1.clear();
	gem2.clear();
	return 0;
    }

    double z_gem1 = Z_gem1;
    double z_gem2 = Z_gem2;
    double z_hycal = Z_hycal;
    double res = Delta;

    vector<GEMClusterStruct> res_gem1;
    vector<GEMClusterStruct> res_gem2;

    int nh = pHHit->size();
    for(int i=0;i<nh;i++)
    {
	double dr = 0.;

	int match_gem1 = 0;
	int match_gem2 = 0;
	float m_index = 0;
	float m_e = 0;

	// search GEM1
	res = Delta;
	int n = gem1.size();
	if(n > 0)
	{
	    double x_hycal = (pHHit->at(i).x) *z_gem1/z_hycal;
	    double y_hycal = (pHHit->at(i).y) *z_gem1/z_hycal;
	    for(int j=0;j<n;j++)
	    {
                dr = TMath::Sqrt((x_hycal - gem1[j].x)*(x_hycal - gem1[j].x)+
		                 (y_hycal - gem1[j].y)*(y_hycal - gem1[j].y));
		if( (dr - res)<=0. )
		{
		    res = dr;
		    m_index = j;
		    m_e = pHHit->at(i).E;
		    match_gem1 = 1;
		    match_gem2 = 0;
		}
	    }
	}

	//continue search GEM2
	res = Delta;
	n = gem2.size();
	if(n>0)
	{
	    double x_hycal = (pHHit->at(i).x) *z_gem2/z_hycal;
	    double y_hycal = (pHHit->at(i).y) *z_gem2/z_hycal;
	    for(int j=0;j<n;j++)
	    {
                dr = TMath::Sqrt((x_hycal - gem2[j].x)*(x_hycal - gem2[j].x)+
		                 (y_hycal - gem2[j].y)*(y_hycal - gem2[j].y));
		if( (dr - res)<=0. )
		{
		    res = dr;
		    m_index = j;
		    m_e = pHHit->at(i).E;
		    match_gem2 = 1;
		    match_gem1 = 0;
		}
	    }
	}
	//after searching
	if( (match_gem1 == 1) && (match_gem2 == 1) ) 
	    cout<<"GEMPhysHandler::HyCalGEMPosMatch(): error...."
	        <<endl;
	else if( match_gem1 == 1 ) 
	{ 
	    gem1[m_index].energy = m_e; 
	    res_gem1.push_back(gem1[m_index]);
	}
	else if( match_gem2 == 1 ) 
	{ 
	    gem2[m_index].energy = m_e; 
	    res_gem2.push_back(gem2[m_index]);
	}
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

#ifndef KEEP_MAX_TWO_CLUSTERS_AFTER_MATCH
    return res_gem1.size() + res_gem2.size();
#endif

#ifdef KEEP_MAX_TWO_CLUSTERS_AFTER_MATCH
    if ( res_gem1.size() + res_gem2.size() > 2)
	return 0;
    else 
	return res_gem1.size() + res_gem2.size();
#endif
}

template<class T> void GEMPhysHandler::EvalMatchMech(T * online_hit)
{
    vector<GEMClusterStruct> gem1, gem2;

    // ep moller ratio before cut using gem
    online_hit->GetClusterHyCalCutMode(gem1, gem2);
    hQuantityOfClustersGEMBeforeMatch -> Fill( gem1.size() + gem2.size() );

    gem1.clear();
    gem2.clear();
    online_hit -> GetClusterHyCalPlusMode(gem1, gem2);

    int nhits_gem1 = gem1.size();
    int nhits_gem2 = gem2.size();
    int nhits_hycal = pHyCalHit->size();
    if(nhits_hycal <=0 )
        return;

    double z_gem1 = Z_gem1;
    double z_gem2 = Z_gem2;
    double z_hycal = Z_hycal;
    double res = Delta; // a larger range, 60mm

    vector<GEMClusterStruct> res_gem1;
    vector<GEMClusterStruct> res_gem2;

    // gem efficiency in production runs
    if( (pHyCalHit->size() == 2) && 
        ( (pHyCalHit->at(0).E + pHyCalHit->at(1).E) > (1000.*beamEnergy - 300.)  ) )
    {
       HyCalMollerElectronQuantity += 2.0;
    }
    else if( (pHyCalHit->size() == 1) && 
             ( pHyCalHit->at(0).E > (1000.*beamEnergy - 300.)  ) )
    {
       //HyCalEpElectronQuantity += 1.0;
    }

    int nh = nhits_hycal;
    for(int i=0;i<nh;i++)
    {
	double dr = 0.;

	int match_gem1 = 0;
	int match_gem2 = 0;
	float m_index = 0;
	float m_e = 0;

	double x_project = 0.;
	double y_project = 0.;

	// search GEM1
	res = Delta;
	int n = gem1.size();
	if(n > 0)
	{
	    double x_hycal = (pHyCalHit->at(i).x)*z_gem1/z_hycal;
	    double y_hycal = (pHyCalHit->at(i).y)*z_gem1/z_hycal;
	    for(int j=0;j<n;j++)
	    {
                dr = TMath::Sqrt((x_hycal - gem1[j].x)*(x_hycal - gem1[j].x)+
		                 (y_hycal - gem1[j].y)*(y_hycal - gem1[j].y));
		if( ( dr - res) <= 0. )
		{
		    res = dr;
		    m_index = j;
		    m_e = pHyCalHit->at(i).E;
		    match_gem1 = 1;
		    match_gem2 = 0;

		    x_project =x_hycal;
		    y_project =y_hycal;
		}
	    }
	}

	//continue search GEM2
	res = Delta;
	n = gem2.size();
	if(n>0)
	{
	    double x_hycal = (pHyCalHit->at(i).x) *z_gem2/z_hycal;
	    double y_hycal = (pHyCalHit->at(i).y) *z_gem2/z_hycal;
	    for(int j=0;j<n;j++)
	    {
                dr = TMath::Sqrt((x_hycal - gem2[j].x)*(x_hycal - gem2[j].x)+
		                 (y_hycal - gem2[j].y)*(y_hycal - gem2[j].y));
		if( (dr - res)<=0. )
		{
		    res = dr;
		    m_index = j;
		    m_e = pHyCalHit->at(i).E;
		    match_gem2 = 1;
		    match_gem1 = 0;

		    x_project = x_hycal;
		    y_project = y_hycal;
		}
	    }
	}
	//after searching
	if( (match_gem1 == 1) && (match_gem2 == 1) ) 
	{
	    cout<<"GEMPhysHandler::HyCalGEMPosMatch(): error...."
	        <<endl;
	    return;
	}
	else if( match_gem1 == 1 ) 
	{ 
	    hXDiffMatch->Fill(gem1[m_index].x - x_project);
	    hYDiffMatch->Fill(gem1[m_index].y - y_project);
	    hhXDiffMatchVsEnergy->Fill(m_e, gem1[m_index].x - x_project);
	    hhYDiffMatchVsEnergy->Fill(m_e, gem1[m_index].y - y_project);

	    gem1[m_index].energy = m_e; 
	    res_gem1.push_back(gem1[m_index]);

	}
	else if( match_gem2 == 1 ) 
	{ 
	    hXDiffMatch->Fill(gem2[m_index].x - x_project);
	    hYDiffMatch->Fill(gem2[m_index].y - y_project);
	    hhXDiffMatchVsEnergy->Fill(m_e, gem2[m_index].x - x_project);
	    hhYDiffMatchVsEnergy->Fill(m_e, gem2[m_index].y - y_project);

	    gem2[m_index].energy = m_e; 
	    res_gem2.push_back(gem2[m_index]);
	    //sectorize
	}
    }

    if(nhits_hycal == 2)
    {
	for(int i=0;i<nhits_hycal;i++)
	{
	    hhHyCalClusterMap2ClusterBeforeMatch->Fill(pHyCalHit->at(i).x, pHyCalHit->at(i).y);
	}
    }

    if( res_gem1.size() + res_gem2.size() == 2 )
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

    if( (res_gem1.size() + res_gem2.size()) > 0)
    {
	if( (nhits_gem1 == 1) && (nhits_gem2==0) &&  (nhits_hycal >=2) )
	{
	    double theta = TMath::Sqrt( (res_gem1[0].x)*(res_gem1[0].x)+(res_gem1[0].y)*(res_gem1[0].y)  )/z_gem1 *180 /PI;
	    hhEnergyVsAngleHyCal2CGEM1C->Fill(theta, gem1[0].energy);
	    hhGEMClusterMapHyCal2GEM1->Fill( res_gem1[0].x, res_gem1[0].y);
	    for(int i=0;i<nhits_hycal;i++)
	    {
		hhHyCalClusterMapHyCal2GEM1->Fill(pHyCalHit->at(i).x, pHyCalHit->at(i).y);
	    }
	}
	else if( (nhits_gem2 == 1) && (nhits_gem1==0) &&  (nhits_hycal >=2) )
	{
	    double theta = TMath::Sqrt( (res_gem2[0].x)*(res_gem2[0].x)+(res_gem2[0].y)*(res_gem2[0].y)  )/z_gem2 *180 /PI;
	    hhEnergyVsAngleHyCal2CGEM1C->Fill(theta, gem2[0].energy);
	    hhGEMClusterMapHyCal2GEM1->Fill( res_gem2[0].x, res_gem2[0].y);
	    for(int i=0;i<nhits_hycal;i++)
	    {
		hhHyCalClusterMapHyCal2GEM1->Fill(pHyCalHit->at(i).x, pHyCalHit->at(i).y);
	    }

	}
    }
    hNbPointsMatch->Fill( res_gem1.size() + res_gem2.size() );
    hQuantityOfClustersGEMAfterMatch->Fill( res_gem1.size() + res_gem2.size() );

    // gem efficiency in production runs
    if( (pHyCalHit->size() == 2) && 
        ( (pHyCalHit->at(0).E + pHyCalHit->at(1).E) > (1000.*beamEnergy - 300.)  ) )
    {
       GEMMollerElectronQuantity += res_gem1.size() + res_gem2.size();
    }
    else if( (pHyCalHit->size() == 1) && 
        ( pHyCalHit->at(0).E > (1000.*beamEnergy - 300.)  ) )
    {
       int nn = res_gem1.size() + res_gem2.size();
       assert( nn <=1 );
       GEMEpElectronQuantity += res_gem1.size() + res_gem2.size();
       HyCalEpElectronQuantity +=1.0;
       
       //sectorize
       double _x_project ;
       double _y_project ;
       float x = pHyCalHit->at(0).x;
       float y = pHyCalHit->at(0).y;
       if(x > 0)
       {
           _x_project = x *z_gem2/z_hycal;
	   _y_project = y *z_gem2/z_hycal;
       }
       else
       {
           _x_project = x *z_gem1/z_hycal;
	   _y_project = y *z_gem1/z_hycal;
       }
       int index = GetSectorIndex(_x_project, _y_project);
       if(index >= 0)
       {
	   //HyCalEpElectronQuantity +=1.0;
           hycal_ep_quantity[index] += 1.0;
           gem_ep_quantity[index] +=  res_gem1.size() + res_gem2.size();
       }

       float r = TMath::Sqrt(x*x+y*y);
       int r_index = (int) (r/10.);
       r_index -= 5;
       if( (r_index >= 0) && (r_index<50) )
       {
           gem_r_sec_ep_quantity[r_index] += res_gem1.size() + res_gem2.size();
	   hycal_r_sec_ep_quantity[r_index] +=1.0;
       }
    }
}


void GEMPhysHandler::GeometryMollerRing(vector<GEMClusterStruct> &gem1, 
                                        vector<GEMClusterStruct>& gem2)
{
    double z_gem1 = Z_gem1;
    double z_gem2 = Z_gem2;
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

double GEMPhysHandler::MollerAngleFromEnergy(double E)
{
    //double E_beam = 2200; // MeV
    double E_beam = beamEnergy * 1000.; //GeV -> MeV
    double E_elec = 0.511; //MeV
    double top = E * E_beam + E_elec*(E - E_beam) - E_elec * E_elec;
    double bottom = TMath::Sqrt( (E_beam * E_beam - E_elec * E_elec) * (E * E - E_elec * E_elec) );
    double cosangle = top / bottom;
    return TMath::ACos(cosangle)*180./3.1415926;
}

double GEMPhysHandler::MollerEnergyFromAngle(double theta)
{
    // Note:
    //      theta in deg.
    //      not in rad
    
    double E_beam = beamEnergy * 1000.;
    double E_elec = 0.511;//MeV
    theta = theta/180. * 3.1415926;
    double _alpha = (E_beam + E_elec) / TMath::Cos(theta) / TMath::Sqrt( E_beam * E_beam - E_elec * E_elec);
    double alpha = _alpha * _alpha;
    return (alpha + 1)/(alpha -1) * E_elec;
}

int GEMPhysHandler::GetSectorIndex(double x, double y)
{
    int index = -1;
    int index_x = -1;
    int index_y = -1;

    if( (x>x_sector[0]) && (x<=x_sector[1]) )
        index_x = 1;
    else if( (x>x_sector[2]) && (x<=x_sector[3]) )
        index_x = 2;
    else if( (x>x_sector[4]) && (x<=x_sector[5]) )
        index_x = 3;
    else if( (x>x_sector[5]) && (x<=x_sector[6]) )
        index_x = 4;
    else if( (x>x_sector[6]) && (x<=x_sector[7]) )
        index_x = 5;
    else if( (x>x_sector[8]) && (x<=x_sector[9]) )
        index_x = 6;
    else if( (x>x_sector[10]) && (x<=x_sector[11]) )
        index_x = 7;
    else if( (x>x_sector[12]) && (x<=x_sector[13]) )
        index_x = 8;

    if( (y>y_sector[0]) && (y<=y_sector[1]) )
        index_y = 1;
    else if( (y>y_sector[2]) && (y<=y_sector[3]) )
        index_y = 2;
    else if( (y>y_sector[4]) && (y<=y_sector[5]) )
        index_y = 3;
    else if( (y>y_sector[5]) && (y<=y_sector[6]) )
        index_y = 4;
    else if( (y>y_sector[7]) && (y<=y_sector[8]) )
        index_y = 5;
    else if( (y>y_sector[8]) && (y<=y_sector[9]) )
        index_y = 6;
    else if( (y>y_sector[10]) && (y<=y_sector[11]) )
        index_y = 7;
    else if( (y>y_sector[12]) && (y<=y_sector[13]) )
        index_y = 8;
    else if( (y>y_sector[14]) && (y<=y_sector[15]) )
        index_y = 9;

    index = (index_y-1) * 8 + index_x;
    if(index > 0 && index <=27 )
        index = index;
    else if(index > 27 && index <=35 )
        index = index -1;
    else if(index > 36 && index <= 72 )
        index = index -2;

    return index -1;
}

void GEMPhysHandler::WriteSectorEff()
{
    outfile<<"<<<>>><<<>>><<<>>><<<>>>"<<endl;
    outfile<<"sector effiency..."<<endl;
    double gem = 0.0;
    double hycal = 0.0;
    /*
    for(int i=0;i<195;i++)
    {
        outfile<<"Moller Efficiency Sector "<<i<<": "
	       <<gem_moller_quantity[i]  <<" / "
	       <<hycal_moller_quantity[i]<<" = "
	       <<gem_moller_quantity[i]/hycal_moller_quantity[i]
	       <<endl;
    }
    outfile<<"......"<<endl;
    */
    for(int i=0;i<70;i++)
    {
       gem += gem_ep_quantity[i];
       hycal += hycal_ep_quantity[i];
       outfile<<"ep Efficiency Sector "<<i<<": "
	       <<gem_ep_quantity[i]  <<" / "
	       <<hycal_ep_quantity[i]<<" = "
	       <<gem_ep_quantity[i]/hycal_ep_quantity[i]
	       <<endl;
    }
    outfile<<"in all "
           <<gem<<" / "
	   <<hycal<<" = "
	   <<gem/hycal
	   <<endl;
    outfile<<"<<<>>><<<>>><<<>>><<<>>>"<<endl;

    outfile<<"<<<>>><<<>>><<<>>><<<>>>"<<endl;
    outfile<<"circular sector effiency..."<<endl;
    double gem1 = 0.0;
    double hycal1 = 0.0;
    for(int i=0;i<50;i++)
    {
       gem1 += gem_r_sec_ep_quantity[i];
       hycal1 += hycal_r_sec_ep_quantity[i];
       outfile<<"ep Efficiency Sector "<<i<<": "
	       <<gem_r_sec_ep_quantity[i]  <<" / "
	       <<hycal_r_sec_ep_quantity[i]<<" = "
	       <<gem_r_sec_ep_quantity[i]/hycal_r_sec_ep_quantity[i]
	       <<endl;
    }
    outfile<<"in all "
           <<gem1<<" / "
	   <<hycal1<<" = "
	   <<gem1/hycal1
	   <<endl;
    outfile<<"<<<>>><<<>>><<<>>><<<>>>"<<endl;

}
