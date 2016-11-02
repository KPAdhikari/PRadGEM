#include "GEMPhysics.h"
#include <iostream>
#include "GEMPedestal.h"
#include "GEMMapping.h"
#include "GEMOnlineHitDecoder.h"
#include <TCanvas.h>
#include "GEMTree.h"
#include "GEMCoord.h"
#include "PRadDataHandler.h"
#include "PRadReconstructor.h"
#include "PRadMoller.h"
#include "PRadEP.h"
#include "HyCalGEMMatch.h"

//#define PROD_GEM_OFFSET

using namespace std;

GEMPhysics::GEMPhysics()
{
    cout<<"GEMPhysics Constructor..."
        <<endl;
    pedestal = new GEMPedestal();
    pedestal -> LoadPedestal();
    mapping = GEMMapping::GetInstance();
    hit_decoder = new GEMOnlineHitDecoder();
    hit_decoder -> SetPedestal(pedestal);

    gem_coord = new GEMCoord();
    gem_coord -> SetHitDecoder(hit_decoder);
    moller_analyzer = new PRadMoller();
    moller_analyzer->SetBeamEnergy(2147);
    gem1_moller_analyzer = new PRadMoller();
    gem1_moller_analyzer->SetBeamEnergy(2147);
    gem2_moller_analyzer = new PRadMoller();
    gem2_moller_analyzer->SetBeamEnergy(2147);
    ep_analyzer = new PRadEP();
    ep_analyzer->SetBeamEnergy(2147);
    gem1_ep_analyzer = new PRadEP();
    gem1_ep_analyzer->SetBeamEnergy(2147);
    gem2_ep_analyzer = new PRadEP();
    gem2_ep_analyzer->SetBeamEnergy(2147);

    match_analyzer = new HyCalGEMMatch();

    InitProductionEfficiency();
}

GEMPhysics::~GEMPhysics()
{
}

void GEMPhysics::SetPRadReconstructor(PRadReconstructor * reconstructor)
{
    reconstruct = reconstructor;
}

void GEMPhysics::SetPRadDataHandler(PRadDataHandler *handler)
{
    pHandler = handler;
}

void GEMPhysics::SetGEMTree(GEMTree *tree)
{
    rst_tree = tree;
}

void GEMPhysics::AccumulateEvent(int evtID, unordered_map<int, vector<int> > & event)
{
    //cout<<"event number from gem: "<<evtID<<endl;
    SetEvtID( evtID );
    hit_decoder -> ProcessEvent(event);
    CharactorizeGEM();
    CharactorizePhysics();
}

void GEMPhysics::CharactorizeGEM()
{
    int n = mapping->GetNbOfDetectors();

    vector<GEMClusterStruct> gem;
    for(int i=0;i<n;i++)
    {
	gem.clear();
	gem_coord->GetClusterGEM(i, gem);
	rst_tree -> PushDetector(i, gem);
    }
    rst_tree -> FillGEMTree();
}

void GEMPhysics::CharactorizePhysics()
{
    hycal_hit = reconstruct->CoarseHyCalReconstruct(pHandler->GetEventCount()-1);
    match_analyzer -> SetGEMCoord(gem_coord);
    match_analyzer -> SetHyCalVector( &hycal_hit);

    //#ifndef PROD_GEM_OFFSET
    CharactorizePlanePhysics();
    //#else
    CharactorizeOverlapPhysics();
    //#endif
}

void GEMPhysics::CharactorizePlanePhysics()
{
    // gem offsets from calibration run
    //gem_coord -> SetGEMOffsetX(0.0674);
    //gem_coord -> SetGEMOffsetY(1.664);
    
    // gem offsets from production run
    //gem_coord -> SetGEMOffsetX(0.3255);
    //gem_coord -> SetGEMOffsetY(0.1924);

    gem_coord -> SetGEMOffsetX(0.);
    gem_coord -> SetGEMOffsetY(0.);
 
    match_analyzer -> Match();
    vector<GEMClusterStruct> gem;
    gem = match_analyzer->GetMatchGEM();

    UpdateProductionEfficiency(gem);

    if(gem.size() == 2){
	if(IsHyCalGEMMatchSucess(gem)){
	    moller_analyzer -> SetEvtID( evt_id);
	    moller_analyzer -> ProcessData(&gem);
	    moller_analyzer -> SetHyCalMatch( match_analyzer->GetMatchHyCal() );
	    rst_tree->PushMoller(moller_analyzer);
	    rst_tree->FillMollerTree();
	}
    }
    else if( gem.size() == 1) {
	ep_analyzer -> SetEvtID(evt_id);
	ep_analyzer -> ProcessData(&gem);
	ep_analyzer -> SetHyCalMatch( match_analyzer->GetMatchHyCal() );
	rst_tree->PushEP(ep_analyzer);
	rst_tree->FillEPTree();
    }
    else{
	// to be implemented...
	// need to select moller or ep 
	// out of the matched more than 3 hits
    }
    moller_analyzer->Reset();
    ep_analyzer->Reset();
    match_analyzer->Reset();
}

void GEMPhysics::CharactorizeOverlapPhysics()
{
    // moller events on single gem
    // to determine gem offsets, need to 
    // set gem offsets to zero
    gem_coord -> SetGEMOffsetX(0.);
    gem_coord -> SetGEMOffsetY(0.);
    match_analyzer->MatchByGEM();
    vector<GEMClusterStruct> gem1;
    vector<GEMClusterStruct> gem2;
    gem1 = match_analyzer -> GetMatchGEM1();
    gem2 = match_analyzer -> GetMatchGEM2();

    if( gem1.size() == 2){
	if(IsHyCalGEMMatchSucess(gem1)){
	    gem1_moller_analyzer -> SetEvtID( evt_id);
	    gem1_moller_analyzer -> ProcessData(&gem1);
	    gem1_moller_analyzer -> SetHyCalMatch( match_analyzer -> GetMatchHyCalGEM1() );
	    rst_tree->PushProdOffset(0, gem1_moller_analyzer);
	    //rst_tree->FillProdOffsetTree();
	}
    }
    else if( gem1.size() == 1){
        gem1_ep_analyzer -> SetChamberID(gem1[0].chamber_id);
	gem1_ep_analyzer -> SetEvtID(evt_id);
	gem1_ep_analyzer -> ProcessData(&gem1);
	gem1_ep_analyzer -> SetHyCalMatch( match_analyzer -> GetMatchHyCalGEM1() );
    }

    if( gem2.size() == 2){
	if(IsHyCalGEMMatchSucess(gem2)){
	    gem2_moller_analyzer -> SetEvtID( evt_id);
	    gem2_moller_analyzer -> ProcessData(&gem2);
	    gem2_moller_analyzer -> SetHyCalMatch( match_analyzer -> GetMatchHyCalGEM2() );
	    rst_tree->PushProdOffset(1, gem2_moller_analyzer);
	    //rst_tree->FillProdOffsetTree();
	}
    }
    else if(gem2.size() == 1){
	gem2_ep_analyzer -> SetChamberID(gem2[0].chamber_id);
	gem2_ep_analyzer -> SetEvtID(evt_id);
	gem2_ep_analyzer -> ProcessData(&gem2);
	gem2_ep_analyzer -> SetHyCalMatch( match_analyzer -> GetMatchHyCalGEM2() );
    }

    if( gem2.size()==2 || gem1.size()==2) 
    {
	rst_tree->FillProdOffsetTree();
    }

    if( gem1.size() ==2 && gem2.size() == 2){
	rst_tree->PushOverlapMollerTree(gem1_moller_analyzer, gem2_moller_analyzer);
	rst_tree->FillOverlapTree();
    }
    else if( gem1.size() == 1 && gem2.size() ==1){
	rst_tree->PushOverlapEpTree(gem1_ep_analyzer, gem2_ep_analyzer);
	rst_tree->FillOverlapTree();
    }

    gem1_moller_analyzer->Reset();
    gem2_moller_analyzer->Reset();
    gem1_ep_analyzer->Reset();
    gem2_ep_analyzer->Reset();
    match_analyzer->Reset();
}

void GEMPhysics::SavePhysResults()
{
    //rst_tree->WriteToDisk();
}

void GEMPhysics::SetEvtID(unsigned int id)
{
    evt_id = id;
}

unsigned int GEMPhysics::GetEvtID()
{
    return evt_id;
}

bool GEMPhysics::IsHyCalGEMMatchSucess(vector<GEMClusterStruct> &gem)
{
    if(gem.size() != 2){
	cout<<" Only can tell moller events."<<endl;
	return true;
    }

    if(gem[0].x == gem[1].x)
	return false;
    else 
	return true;
}

void GEMPhysics::InitProductionEfficiency()
{
    eff_log.open("production_efficiency.log",iostream::out);
    if(!eff_log.is_open()){
        cout<<"cannot open eff log file."<<endl;
	exit(-1);
    }

    for(int i=0;i<21;i++){
	x_hycal_ep_quantity[i]=0.;
	x_gem_ep_quantity[i]=0.;
	y_hycal_ep_quantity[i]=0.;
	y_gem_ep_quantity[i]=0.;
    }

    gem_ep_quantity = 0.;
    hycal_ep_quantity = 0.;
}

static int get_sector_index(float val)
{
    int index = -1;

    if(val>=-600. && val<-552.86)
	index = 0;
    else if(val>=-532.86 && val <-495.72)
	index = 1;
    else if(val>=-475.72 && val <-438.58)
	index = 2;
    else if(val>=-418.58 && val <-381.44)
	index = 3;
    else if(val>=-361.44 && val <-324.30)
	index = 4;
    else if(val>=-304.30 && val <-267.16)
	index = 5;
    else if(val>=-247.16 && val <-210.02)
	index = 6;
    else if(val>=-190.02 && val <-152.88)
	index = 7;
    else if(val>=-132.88 && val <-95.74)
	index = 8;
    else if(val>=-75.74 && val <-38.60)
	index = 9;
    else if(val>=-18.6 && val <18.6)
	index = 10;
    else if(val>=38.6 && val <75.74)
	index = 11;
    else if(val>=95.74 && val <132.88)
	index = 12;
    else if(val>=152.88 && val <190.02)
	index = 13;
     else if(val>=210.02 && val <247.16)
	index = 14;
    else if(val>=267.16 && val <304.30)
	index = 15;
     else if(val>=324.30 && val <361.44)
	index = 16;
    else if(val>=381.44 && val <418.58)
	index = 17;
    else if(val>=438.58 && val <475.72)
	index = 18;
     else if(val>=495.72 && val <532.86)
	index = 19;
    else if(val>=552.86 && val <600)
	index = 20;

    return index;
}

void GEMPhysics::UpdateProductionEfficiency(vector<GEMClusterStruct> &gem)
{
    int n_hycal = hycal_hit.size();
    float tot_e = 0.;
    int x_index = -1;
    int y_index = -1;

    if(n_hycal!=1)
	return;

    tot_e += hycal_hit[0].E;
    if(tot_e>2500. || tot_e<1900.)
	return;

    int n_gem = gem.size();
    if(n_gem > 1){
	cout<<"something very wrong happened..."<<endl;
    }

    x_index = get_sector_index(hycal_hit[0].x *5260./5800.);
    y_index = get_sector_index(hycal_hit[0].y *5260./5800.);

    hycal_ep_quantity+=1.0;
    gem_ep_quantity += n_gem;
    if(x_index != -1){
	x_hycal_ep_quantity[ x_index ]+=1.0;
	x_gem_ep_quantity[ x_index ]+=n_gem;
    }
    if(y_index!=-1){
	y_hycal_ep_quantity[ y_index ]+=1.0;
	y_gem_ep_quantity[ y_index ]+=n_gem;
    }

    if((int)hycal_ep_quantity%1000 == 0 && hycal_ep_quantity>100)
    {
	eff_log<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
	eff_log<<"overall efficiency: "<<gem_ep_quantity<<"/"<<hycal_ep_quantity<<"="
	    <<gem_ep_quantity/hycal_ep_quantity<<endl;
	for(int i=0;i<21;i++){
	    eff_log<<"x efficiency: "<<x_gem_ep_quantity[i]<<"/"<<x_hycal_ep_quantity[i]<<"="
		<<x_gem_ep_quantity[i]/x_hycal_ep_quantity[i]<<endl;
	}
	for(int i=0;i<21;i++){
	    eff_log<<"y efficiency "<<y_gem_ep_quantity[i]<<"/"<<y_hycal_ep_quantity[i]<<"="
		<<y_gem_ep_quantity[i]/y_hycal_ep_quantity[i]<<endl;
	}
    }
}
