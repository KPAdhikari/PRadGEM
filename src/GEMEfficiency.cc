//============================================================================//
// Analyze GEM Efficiency From Hycal Calibration and Gain Equalizing          //
//                                                                            //
// Xinzhan Bai                                                                //
// 08/14/2016                                                                 //
//============================================================================//
#include "GEMEfficiency.h"
#include "GEMConfigure.h"
#include "GEMTree.h"
#include <iostream>
#include "GEMCoord.h"
#include "GEMOnlineHitDecoder.h"
#include "GEMPedestal.h"
#include "GEMMapping.h"
#include "GEMDataStruct.h"
#include "HyCalGEMMatch.h"
#include "PRadDataHandler.h"
#include "PRadReconstructor.h"

//#define POSITION_MATCH_HYCAL

#define SAVE_CALIBRATION_DATA_TREE

using namespace std;

GEMEfficiency::GEMEfficiency()
{
    mapping = GEMMapping::GetInstance();
    pedestal = new GEMPedestal();
    pedestal -> LoadPedestal();

    hit_decoder = new GEMOnlineHitDecoder();
    hit_decoder -> SetPedestal(pedestal);
    gem_coord = new GEMCoord();
    gem_coord -> SetHitDecoder(hit_decoder);

    match_analyzer = new HyCalGEMMatch();
    match_analyzer -> SetMatchCriteria(60.);

    period_trigger = 0.;
    period_gem = 0.;
    n_trigger = 0.;
    n_gem = 0.;
#ifdef CHECK_PAIR_COMPTON
    compton_event = 0.;
    total_event = 0.;
    pair_event = 0.;
#endif
}

GEMEfficiency::~GEMEfficiency()
{   
    cout<<"GEMEfficiency destructor."<<endl;
    log_file<<"overall trigger: "
            <<n_trigger
	    <<", overall gem: "
	    <<n_gem
	    <<", overall eff."
	    <<n_gem/n_trigger
	    <<endl;
#ifdef CHECK_PAIR_COMPTON
    cout<<"total event:"
	<<total_event
	<<", compton event:"
	<<compton_event
	<<", pair event:"
	<<pair_event
	<<endl;
    log_file<<"//////////////////////////"
	<<endl;
    log_file<<"total event:"
	<<total_event
	<<", compton event:"
	<<compton_event
	<<", pair event:"
	<<pair_event
	<<endl;
#endif
}

void GEMEfficiency::AccumulateEvent(unsigned int evtID, unordered_map<int, vector<int> > & event,
	unordered_multimap<string, double> & tdc_channel_value)
{
    SetEventNumber(evtID);

#ifdef SAVE_CALIBRATION_DATA_TREE
    hit_decoder->ProcessEvent(event);

    vector<GEMClusterStruct> gem1, gem2;
    //gem_coord -> GetClusterGEMPlusMode(0, gem1);
    //gem_coord -> GetClusterGEMPlusMode(1, gem2);
    gem_coord -> GetClusterGEM(0, gem1);
    gem_coord -> GetClusterGEM(1, gem2);


    hycal_hit = reconstructor -> CoarseHyCalReconstruct(pHandler->GetEventCount()-1);
    gem_tree->PushCalibrationData(gem1, gem2, tdc_channel_value, &hycal_hit);
    gem_tree->FillGEMTree();

#else
    if( TDCCut(tdc_channel_value) )
    {
	hit_decoder->ProcessEvent(event);

	SaveGEM();
	AlignGEM();
	if( passHyCalGEMPosMatch() )
	{
	    IncreasePeriod(1.0);
	}
	else {
	    IncreasePeriod(0.);
	}
    }
#endif
}

void GEMEfficiency::InitTDCCutVar()
{
    hycal_group_q = configure->TDC_Quan;
    for(int i=0;i<hycal_group_q;i++){
	hycal_group[i] = configure->TDC[i];
    }
    scin_mode = configure->TDC_Channel;

    if(scin_mode == "S1") {
	scin_group[0] = "S1";
	scin_q = 1;
    }
    else if(scin_mode == "S2"){
	scin_group[0] = "S2";
	scin_q = 1;
    }
    else if(scin_mode == "S1orS2") {
	scin_group[0] = "S1";
	scin_group[1] = "S2";
	scin_q = 2;
    }
    else if(scin_mode == "S1andS2") {
	scin_group[0] = "S1";
	scin_group[1] = "S2";
	scin_q = 3;
    }
    else {
	cout<<" ERROR: TDC Cut: unkonw scin. cut mode."<<endl;
    }

    hycal_start = configure->Hycal_Timing_Cut_Start;
    hycal_end = configure->Hycal_Timing_Cut_End;
    scin_start = configure->TDC_Start;
    scin_end = configure->TDC_End;

    log_file<<"--- calibration efficiency cut setting ---"<<endl;
    for(int i=0;i<hycal_group_q;i++){
	log_file<<"hycal group: "<<hycal_group[i]<<endl;
    }
    log_file<<"scin cut mode: "<<scin_mode<<endl;
    log_file<<"hycal tdc start: "<<hycal_start<<" end: "<<hycal_end<<endl;
    log_file<<"scin  tdc start: "<<scin_start<<" end: "<<scin_end<<endl;
    log_file<<"--- calibration efficiency cut setting ---"<<endl;
}

void GEMEfficiency::SaveGEM()
{
    int n = mapping->GetNbOfDetectors();

    vector<GEMClusterStruct> gem;
    for(int i=0;i<n;i++){
	gem.clear();
	gem_coord -> GetClusterGEM(i, gem);
	gem_tree -> PushDetector(i, gem);
    }
    gem_tree -> FillGEMTree();
}

bool GEMEfficiency::TDCCut(unordered_multimap<string, double> & tdc_map)
{
    //------------------------------------
    //  This is to be improved...
    //      for now, the method is:
    //
    //           scin. has signal whose timing
    //           is in some range (min, max)
    //
    //           hycal has signal whose timing
    //           is in some range (min, max)
    //
    //           and use AND for the above two
    //
    //       Better to be:
    //           T_hycal - T_scin is within
    //           some range (min, max)
    //
    //           to be implemented...
    //------------------------------------

    bool scin_b = ScinTimingCut(tdc_map);
    bool hycal_b = HyCalTimingCut(tdc_map);
    return scin_b && hycal_b;
}

bool GEMEfficiency::ScinTimingCut(unordered_multimap<string, double> &tdc_map)
{
    bool scin_b = false;
    vector<float> scin_val[2];
    int nit = 0;
    if( scin_q <= 3  && scin_q > 0){
	if( scin_q == 3)
	    nit = 2;
	else
	    nit = scin_q;
    }
    else
	cout<<" Scin. timing cut: unknown mode."<<endl;

    for(int i=0;i<nit;i++){
	auto range = tdc_map.equal_range( (string)scin_group[i] );
	for(auto it=range.first;it!=range.second;++it)
	    scin_val[i].push_back( (*it).second );
    }

    if( scin_q < 3) // or mode
    {
	for(int i=0;i<scin_q;i++){
	    for( auto &ii: scin_val[i])
		if( ii>scin_start && ii<scin_end)
		    scin_b = true;
	}
    }
    else if( scin_q == 3) // and mode
    {
	bool b1 = false;
	bool b2 = false;
	for( auto &i: scin_val[0])
	    if( i>scin_start && i<scin_end)
		b1 = true;
	for( auto &j: scin_val[1])
	    if( j>scin_start && j<scin_end)
		b2 = true;
	scin_b = b1 && b2;
    }
    else
	cout<<" scin. unkown cut mode."<<endl;
    return scin_b;
}

bool GEMEfficiency::HyCalTimingCut(unordered_multimap<string, double> &tdc_map)
{
    bool hycal_b = false;
    vector<float> hycal_val[4];
    for(int i=0;i<hycal_group_q;i++){
	auto range = tdc_map.equal_range( (string)hycal_group[i] );
	for(auto it=range.first;it!=range.second;++it)
	    hycal_val[i].push_back( (*it).second );
    }

    // hycal OR mode
    for(int i=0;i<hycal_group_q;i++){
	for(auto &ii: hycal_val[i] )
	    if( ii>hycal_start && ii<hycal_end)
		hycal_b = true;
    }
    return hycal_b;
}

void GEMEfficiency::IncreasePeriod(int i)
{
    if( period_trigger != TRACK_PERIOD){
	period_trigger += 1.0;
	if( i != 0)
	    period_gem += 1.0;
    }
    else {
	log_file<<" period trigger: "
	    <<period_trigger
	    <<", period gem: "
	    <<period_gem
	    <<endl;
	log_file<<" period efficiency: "
	    << period_gem / period_trigger
	    <<endl;
	log_file<<" accumulated trigger: "
	    <<n_trigger
	    <<", accumulated gem: "
	    <<n_gem
	    <<endl;
	log_file<<" accumulated eff. : "
	    <<n_gem / n_trigger
	    <<endl;
	period_trigger = 1.0;
	if( i != 0)
	    period_gem = 1.0;
	else
	    period_gem = 0.;
    }

    IncreaseOverall(i);
}

void GEMEfficiency::IncreaseOverall(int i)
{
    n_trigger += 1.0;
    if( i!= 0)
	n_gem += 1.0;
}

void GEMEfficiency::SetGEMTree(GEMTree* tree)
{
    gem_tree = tree;
}

void GEMEfficiency::SetEventNumber(unsigned int num)
{
    event_number = num;
}

unsigned int GEMEfficiency::GetEventNumber()
{
    return event_number;
}

void GEMEfficiency::SetGEMConfigure(GEMConfigure *config)
{
    configure = config;
    InitTDCCutVar();
    InitLogFile();
}

void GEMEfficiency::SetGEMCoord(GEMCoord* coord)
{
    gem_coord = coord;
}

void GEMEfficiency::SetPRadReconstructor(PRadReconstructor * re)
{
    reconstructor = re;
}

void GEMEfficiency::SetPRadDataHandler(PRadDataHandler * h)
{
    pHandler = h;
}

bool GEMEfficiency::passHyCalGEMPosMatch()
{
    vector<GEMClusterStruct> gem;
#ifdef POSITION_MATCH_HYCAL
    hycal_hit = reconstructor -> CoarseHyCalReconstruct(pHandler->GetEventCount()-1);
    match_analyzer -> SetGEMCoord(gem_coord);
    match_analyzer -> SetHyCalVector( &hycal_hit);
    match_analyzer -> Match();

    gem = match_analyzer->GetMatchGEM();

    if(gem.size() != 0)
	return true;
    else 
	return false;
#else
    int size = 0;
    for(int i=0;i<2;i++){
	gem_coord->GetClusterGEM(i, gem);
	size += gem.size();
    }
    if( size != 0)
	return true;
    else 
	return false;
#endif
}

void GEMEfficiency::InitLogFile()
{
    string txt = configure->phys_results_path;
    txt = txt + string(".log");
    log_file.open(txt.c_str(), ios::out);
    log_file.close();
    log_file.open(txt.c_str(), ios::in|ios::out);
    if(!log_file.is_open() )
	cout<<"Error: open log file failed."
	    <<endl;
}

void GEMEfficiency::AlignGEM()
{
    vector<GEMClusterStruct> gem1;
    vector<GEMClusterStruct> gem2;

#ifdef CHECK_PAIR_COMPTON
    int nth = 0;
    gem_coord->GetClusterGEM(nth, gem1);
    nth = 1;
    gem_coord->GetClusterGEM(nth, gem2);
    total_event +=1.0;
    if( gem1.size() == 1 && gem2.size() == 1)
	compton_event += 1.0;
    else if( gem1.size() ==2 || gem2.size() == 2)
	pair_event += 1.0;
    gem1.clear();
    gem2.clear();
#endif

    hycal_hit = reconstructor -> CoarseHyCalReconstruct(pHandler->GetEventCount()-1);

    gem_coord -> SetGEMOffsetX(0.);
    gem_coord -> SetGEMOffsetY(0.);
    match_analyzer -> SetGEMCoord(gem_coord);
    match_analyzer -> SetHyCalVector( &hycal_hit);
    match_analyzer -> MatchByGEM(); 
    // MatchByGEM(): criteria is hard-coded to 10mm, refer to code

    gem1 = match_analyzer -> GetMatchGEM1();
    gem2 = match_analyzer -> GetMatchGEM2();
    //-------------------------------------------------
    //  MatchByGEM() matching mechanism will only keep 
    //  one closest
    //  cluster on gem. because in calibration run, hycal
    //  only have one cluster, (hycal cannot distinguish
    //  e+ from e- because its position resolution ~ 2cm, 
    //  e+ and e- are very close)

    //--------------------------------------------------
    //  only use one electron event, is this reasonable?
    //
    //  given that e+ and e- should be very close, maybe 
    //  should use all of them? how to match clusters on
    //  two gems?
    //
    //  step 1), match hits between two chambers? wrong.
    //      Must match gem to hycal, otherwise meaningless.
    //      because clusters reconstructed by gem is not 
    //      reliable (minimum ionization).
    //--------------------------------------------------
    if(gem1.size() != 1 || gem2.size() != 1)
	return;

    if( gem1[0].x>-22 && gem1[0].x<22 && gem2[0].x>-22 && gem2[0].x<22)
    {
	double x_offset = -(gem1[0].x - gem2[0].x);
	double y_offset = -(gem1[0].y - gem2[0].y);
	gem_tree->PushCaliOffset(x_offset, y_offset);
	gem_tree->FillCaliOffsetTree();
    }
}
