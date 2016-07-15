#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include "PRadGEMTree.h"
#include <iostream>

PRadGEMTree::PRadGEMTree()
{
    std::cout<<"PRadGEMTree Constructor..."<<std::endl;
    tree = new TTree("PRadGEMTree", "PRadGEMTree");
    tree -> Branch("evt_id", &evt_id, "evt_id/l");
   

    GEMClusterStruct g1; 
    a_gem1 = new TClonesArray("GEMClusterStruct");
    a_gem2 = new TClonesArray("GEMClusterStruct");
    a_hycal = new TClonesArray("HyCalClusterStruct");

    tree -> Branch("gem1", &a_gem1);
    tree -> Branch("gem2", &a_gem2);
    tree -> Branch("hycal", &a_hycal);
}

PRadGEMTree::~PRadGEMTree()
{}

void PRadGEMTree::Fill(const unsigned long &index, const std::vector<GEMClusterStruct> &gem1 , const std::vector<GEMClusterStruct> &gem2, const std::vector<HyCalHit> &hycal)
{
    a_gem1->Clear();
    a_gem2->Clear();
    a_hycal->Clear();
    evt_id = index;

    int i = 0;

    for(i=0; i< gem1.size();i++)
    {
        new((*a_gem1)[i]) GEMClusterStruct(gem1[i]);
    }
    for(i=0;i< gem2.size();i++)
    {
        new((*a_gem2)[i]) GEMClusterStruct(gem2[i]);
    }
    for(i=0;i< hycal.size(); i++)
    {
       new((*a_hycal)[i]) HyCalClusterStruct(hycal[i].x, hycal[i].y, hycal[i].E); 
    }

    tree->Fill();
}

void PRadGEMTree::Save()
{
    TFile *file = new TFile("PRadGEMTree.root","RECREATE");
    tree->Write();
    file->Close(); 
}
