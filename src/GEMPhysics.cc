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
    if(gem.size() == 2){
	if(IsHyCalGEMMatchSucess(gem)){
	    moller_analyzer -> SetEvtID( evt_id);
	    moller_analyzer -> ProcessData(&gem);
	    rst_tree->PushMoller(moller_analyzer);
	    rst_tree->FillMollerTree();
	}
    }
    else if( gem.size() == 1) {
	ep_analyzer -> SetEvtID(evt_id);
	ep_analyzer -> ProcessData(&gem);
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
	    rst_tree->PushProdOffset(0, gem1_moller_analyzer);
	    //rst_tree->FillProdOffsetTree();
	}
    }
    else if( gem1.size() == 1){
        gem1_ep_analyzer -> SetChamberID(gem1[0].chamber_id);
	gem1_ep_analyzer -> SetEvtID(evt_id);
	gem1_ep_analyzer -> ProcessData(&gem1);
    }

    if( gem2.size() == 2){
	if(IsHyCalGEMMatchSucess(gem2)){
	    gem2_moller_analyzer -> SetEvtID( evt_id);
	    gem2_moller_analyzer -> ProcessData(&gem2);
	    rst_tree->PushProdOffset(1, gem2_moller_analyzer);
	    //rst_tree->FillProdOffsetTree();
	}
    }
    else if(gem2.size() == 1){
        gem2_ep_analyzer -> SetChamberID(gem2[0].chamber_id);
	gem2_ep_analyzer -> SetEvtID(evt_id);
	gem2_ep_analyzer -> ProcessData(&gem2);
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
