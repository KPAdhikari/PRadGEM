#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include "PRadGEMTree.h"
#include <iostream>

PRadGEMTree::PRadGEMTree()
{
   // variable init
    Clear();

    // tree
    std::cout<<"PRadGEMTree Constructor..."<<std::endl;
    file = new TFile("./results/PRadGEMTree.root","RECREATE");
    tree = new TTree("PRadGEMTree", "PRadGEMTree");
    tree -> SetDirectory(file);
    tree -> Branch("evt_id", &evt_id, "evt_id/l");

#ifdef CLONE_ARRAY
    //GEMClusterStruct g1; 
    a_gem1 = new TClonesArray("GEMClusterStruct");
    a_gem2 = new TClonesArray("GEMClusterStruct");
    a_hycal = new TClonesArray("HyCalClusterStruct");
    tree -> Branch("gem1", &a_gem1);
    tree -> Branch("gem2", &a_gem2);
    tree -> Branch("hycal", &a_hycal);
#else
    tree -> Branch("gem1_nhit", &gem1_nhit, "gem1_nhit/I");
    tree -> Branch("gem1_x", &gem1_x, "gem1_x[gem1_nhit]/F");
    tree -> Branch("gem1_y", &gem1_y, "gem1_y[gem1_nhit]/F");
    tree -> Branch("gem1_x_charge", &gem1_x_charge, "gem1_x_charge[gem1_nhit]/F");
    tree -> Branch("gem1_y_charge", &gem1_y_charge, "gem1_y_charge[gem1_nhit]/F");
    tree -> Branch("gem1_energy", &gem1_energy, "gem1_energy[gem1_nhit]/F");
    tree -> Branch("gem1_z", &gem1_z, "gem1_z[gem1_nhit]/F");
    tree -> Branch("gem1_x_size", &gem1_x_size, "gem1_x_size[gem1_nhit]/I");
    tree -> Branch("gem1_y_size", &gem1_y_size, "gem1_y_size[gem1_nhit]/I");

    tree -> Branch("gem2_nhit", &gem2_nhit, "gem2_nhit/I");
    tree -> Branch("gem2_x", &gem2_x, "gem2_x[gem2_nhit]/F");
    tree -> Branch("gem2_y", &gem2_y, "gem2_y[gem2_nhit]/F");
    tree -> Branch("gem2_x_charge", &gem2_x_charge, "gem2_x_charge[gem2_nhit]/F");
    tree -> Branch("gem2_y_charge", &gem2_y_charge, "gem2_y_charge[gem2_nhit]/F");
    tree -> Branch("gem2_energy", &gem2_energy, "gem2_energy[gem2_nhit]/F");
    tree -> Branch("gem2_z", &gem2_z, "gem2_z[gem2_nhit]/F");
    tree -> Branch("gem2_x_size", &gem2_x_size, "gem2_x_size[gem2_nhit]/I");
    tree -> Branch("gem2_y_size", &gem2_y_size, "gem2_y_size[gem2_nhit]/I");
    
    tree -> Branch("hycal_nhit", &hycal_nhit, "hycal_nhit/I");
    tree -> Branch("hycal_x", &hycal_x, "hycal_x[hycal_nhit]/F");
    tree -> Branch("hycal_y", &hycal_y, "hycal_y[hycal_nhit]/F");
    tree -> Branch("hycal_energy", &hycal_energy, "hycal_energy[hycal_nhit]/F");
#endif

    // independent physics object
    tree -> Branch("angular_resolution_from_moller1", &angular_resolution_from_moller1, "angular_resolution_from_moller1/F");
    tree -> Branch("angular_resolution_from_moller2", &angular_resolution_from_moller2, "angular_resolution_from_moller2/F");
    tree -> Branch("angle_moller1", &angle_moller1, "angle_moller1/F");
    tree -> Branch("angle_moller2", &angle_moller2, "angle_moller2/F");
    tree -> Branch("symm_dx", &symm_dx, "symm_dx/F");
    tree -> Branch("symm_dy", &symm_dy, "symm_dy/F");
    tree -> Branch("symm_coplanarity", &symm_coplanarity, "symm_coplanarity/F");
}

PRadGEMTree::~PRadGEMTree()
{}

void PRadGEMTree::Clear()
{    
    angular_resolution_from_moller1 = 360.;
    angular_resolution_from_moller2 = 360.;
    angle_moller1 = 360.;
    angle_moller2 = 360.;
    symm_dx = 1200.;
    symm_dy = 1200.;
    symm_coplanarity = 180.;
}


void PRadGEMTree::Fill(const unsigned long &index, 
                       std::vector<GEMClusterStruct> &gem1 , 
		       std::vector<GEMClusterStruct> &gem2, 
		       std::vector<HyCalHit> &hycal)
{
    evt_id = index;

#ifdef CLONE_ARRAY
    a_gem1 -> ExpandCreate(gem1.size());
    a_gem2 -> ExpandCreate(gem2.size());
    a_hycal -> ExpandCreate(hycal.size());
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
    a_gem1->Clear("C");
    a_gem2->Clear("C");
    a_hycal->Clear("C");
#else
    int i = 0;
    gem1_nhit = gem1.size();
    for(i=0; i< gem1.size();i++)
    {
       gem1_x[i] = gem1[i].x; 
       gem1_y[i] = gem1[i].y; 
       gem1_x_charge[i] = gem1[i].x_charge; 
       gem1_y_charge[i] = gem1[i].y_charge; 
       gem1_energy[i] = gem1[i].energy; 
       gem1_z[i] = gem1[i].z; 
       gem1_x_size[i] = gem1[i].x_size; 
       gem1_y_size[i] = gem1[i].y_size; 
 
    }

    gem2_nhit = gem2.size();
    for(i=0;i< gem2.size();i++)
    {
       gem2_x[i] = gem2[i].x; 
       gem2_y[i] = gem2[i].y; 
       gem2_x_charge[i] = gem2[i].x_charge; 
       gem2_y_charge[i] = gem2[i].y_charge; 
       gem2_energy[i] = gem2[i].energy; 
       gem2_z[i] = gem2[i].z; 
       gem2_x_size[i] = gem2[i].x_size; 
       gem2_y_size[i] = gem2[i].y_size; 
 
    }

    hycal_nhit = hycal.size();
    for(i=0;i< hycal.size(); i++)
    {
        hycal_x[i] = hycal[i].x;
	hycal_y[i] = hycal[i].y;
	hycal_energy[i] = hycal[i].E;
    }
    tree->Fill();
#endif
 
}

void PRadGEMTree::Save()
{
    //tree->GetCurrentFile()->Write();
    file->Write();
    //file->Save();
}
