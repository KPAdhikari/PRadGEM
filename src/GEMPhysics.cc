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

#ifndef PROD_GEM_OFFSET
    // gem offsets from calibration run
    //gem_coord -> SetGEMOffsetX(0.0674);
    //gem_coord -> SetGEMOffsetY(1.664);
    gem_coord -> SetGEMOffsetX(0.);
    gem_coord -> SetGEMOffsetY(0.);
    match_analyzer -> Match();
    vector<GEMClusterStruct> gem;
    gem = match_analyzer->GetMatchGEM();
    if(gem.size() == 2){
	moller_analyzer -> ProcessData(&gem);
	rst_tree->PushMoller(moller_analyzer);
	rst_tree->FillMollerTree();
    }
    else if( gem.size() == 1) {
	ep_analyzer -> ProcessData(&gem);
	rst_tree->PushEP(ep_analyzer);
	rst_tree->FillEPTree();
    }
    else{
	// to be implemented...
	// need to select moller or ep 
	// out of the matched more than 3 hits
    }
#else
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
	gem1_moller_analyzer -> ProcessData(&gem1);
	rst_tree->PushProdOffset(0, gem1_moller_analyzer);
	//rst_tree->FillProdOffsetTree();
    }
    if( gem2.size() == 2){
	gem2_moller_analyzer -> ProcessData(&gem2);
	rst_tree->PushProdOffset(1, gem2_moller_analyzer);
	//rst_tree->FillProdOffsetTree();
    }
    if( gem2.size()==2 || gem1.size()==2 )
	rst_tree->FillProdOffsetTree();
#endif
}

void GEMPhysics::SavePhysResults()
{
    //rst_tree->WriteToDisk();
}
