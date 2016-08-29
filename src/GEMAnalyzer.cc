//==============================================================================//
// Interface class, no real work done here                                      //
// however, worker class must be passed here.                                   //
//        e.g.  class GEMEventAnalyzer analyze GEM data.                        //
//                   all GEM related analysis class must be registered to       //
//                   GEMEventAnalyzer                                           //
//              class HyCalEvnetAnalyzer analyze hycal data.                    //
//                   all hycal related analysis class must be registered to     //
//                   HyCalEventAnalyzer                                         //
//              class EpicsEventAnalyzer analyzer epics data.                   //
//                   all epics related analysis class must be registered to     //
//                   EpicsEventAnalyzer                                         //
//              class for other detectors ....                                  //
//                                                                              //
// Xinzhan Bai                                                                  //
// 07/22/2016                                                                   //
//==============================================================================//
#include "GEMAnalyzer.h"
#include "GEMConfigure.h"
#include "GEMEvioParser.h"
#include "PRadBenchMark.h"
#include "GEMDataHandler.h"
#include "GEMEventAnalyzer.h"
#include "GEMPedestal.h"
#include "GEMPhysics.h"
#include "PRadDataHandler.h"
#include "PRadReconstructor.h"
#include "EpicsEventAnalyzer.h"
#include "EpicsPhysics.h"
#include "GEMTree.h"
#include "TDCEventAnalyzer.h"
#include "GEMEfficiency.h"
#include "EventUpdater.h"

static string trim(const string &str, const string &w = " \t\n\r")
{

    const auto strBegin = str.find_first_not_of(w);
    if (strBegin == string::npos)
	return ""; // no content

    const auto strEnd = str.find_last_not_of(w);

    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

GEMAnalyzer::GEMAnalyzer()
{
    res_tree = new GEMTree();

    configure = new GEMConfigure();
    configure->LoadConfigure();

    parser = new GEMEvioParser();
    handler = new GEMDataHandler();
    gem_event_analyzer = new GEMEventAnalyzer();
    pedestal = new GEMPedestal();
    physics = new GEMPhysics();
    epics_event_analyzer = new EpicsEventAnalyzer();
    tdc_event_analyzer = new TDCEventAnalyzer();
    gem_efficiency = new GEMEfficiency();
    gem_efficiency ->SetGEMConfigure(configure);
    gem_efficiency ->SetGEMTree(res_tree);
    event_updater = new EventUpdater();
    epics_physics = new EpicsPhysics();
    GetRunType();
    RegisterClass();
    InitPRadDataHandler();
}

GEMAnalyzer::~GEMAnalyzer()
{
    delete gem_efficiency;
}

void GEMAnalyzer::RegisterClass()
{
    // gem
    physics -> SetGEMTree(res_tree);
    gem_efficiency->SetGEMTree(res_tree);
    event_updater->SetGEMEventAnalyzer(gem_event_analyzer);
    event_updater->SetGEMEfficiencyAnalyzer(gem_efficiency);
    event_updater->SetGEMPhysics(physics);
    event_updater->SetGEMPedestal(pedestal);
    //TDC
    event_updater->SetTDCEventAnalyzer(tdc_event_analyzer);
    // hycal
    // epics
    epics_physics -> SetGEMTree(res_tree);
    event_updater -> SetEpicsEventAnalyzer(epics_event_analyzer);
    event_updater -> SetEpicsPhysics(epics_physics);
    //others ...

    // data handler
    SetDataHandler();
}

void GEMAnalyzer::SetDataHandler()
{
    parser->SetDataHandler(handler);
    parser->SetEventUpdater(event_updater);

    // gem
    handler -> SetGEMEventAnalyzer(gem_event_analyzer);
    // epics
    handler -> SetEpicsEventAnalyzer(epics_event_analyzer);
    // tdc
    handler -> SetTDCEventAnalyzer(tdc_event_analyzer);
    // ...
}

void GEMAnalyzer::InitPRadDataHandler()
{
    // ...
    pHandler = new PRadDataHandler();
    pHandler->ReadTDCList("/home/xbai/w/pRad/source/PRadDecoder/config/tdc_group_list.txt");
    pHandler->ReadChannelList("/home/xbai/w/pRad/source/PRadDecoder/config/module_list.txt");
    pHandler->BuildChannelMap();
    pHandler->ReadPedestalFile("/home/xbai/w/pRad/source/PRadDecoder/config/pedestal.dat");
    pHandler->ReadCalibrationFile("/home/xbai/w/pRad/source/PRadDecoder/config/calibration.txt");

    reconstruct = new PRadReconstructor();
    reconstruct->SetHandler(pHandler);

    // GEMEvioParser
    parser->SetPRadDataHandler(pHandler);
    physics->SetPRadReconstructor(reconstruct);
    physics->SetPRadDataHandler(pHandler);
    gem_efficiency -> SetPRadReconstructor(reconstruct);
    gem_efficiency -> SetPRadDataHandler(pHandler);
}

void GEMAnalyzer::ProcessFiles()
{
    cout<<"## number of files to be processed:  "
	<<configure->nFile <<" ## "
	<<endl;

    int nF = configure->nFile;

    for(int i=0;i<nF;++i)
    {
	file = configure->fileList[i];

	if(gem_physics != 1){
	    // pedestal run or raw data run
	    // these runs do not need hycal information
	    // no need to initialize hycal
	    
	    PRadBenchMark timer;
	    Process();
	    cout<<"used time: "<<timer.GetElapsedTime()<<" ms"<<endl;
	}
	else {
	    // physics run

	    // initialize by data
	    if(file.find("evio.0")!=string::npos){
		pHandler->InitializeByData(configure->fileList[0].c_str());
	    }
	    else {
		PRadBenchMark timer;
		Process();
		cout<<"used time: "<<timer.GetElapsedTime()<<" ms"<<endl;
	    }
	}
    }

    ProcessResults();
}

void GEMAnalyzer::Process()
{
    parser->ParseFile(file);
}

void GEMAnalyzer::ProcessResults()
{
    if(gem_pedestal)
	pedestal -> SavePedestalFile();
    else if(gem_physics)
	;//physics -> SavePhysResults();

    //epics_physics -> SaveEpicsResults();
    res_tree->WriteToDisk();
}

void GEMAnalyzer::GetRunType()
{
    gem_physics = 0;
    gem_pedestal = 0;
    gem_raw = 0;

    string runtype = trim(configure->GetRunType());
    cout<<"### RUN TYPE ### ";
    cout<<configure->GetRunType()<<endl;
    gem_event_analyzer->SetRunType(runtype);
    event_updater -> SetRunType(runtype);
    if(runtype == "PHYSICS")
	gem_physics = 1;
    else if(runtype == "PEDESTAL")
	gem_pedestal = 1;
    else if(runtype == "RAWDATA")
	gem_raw = 1;
    else
    {
	cout<<" GEMAnalyzer: undefined run type..."
	    <<endl;
	return;
    }
}
