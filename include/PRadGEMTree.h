#ifndef _PRAD_GEM_TREE_H
#define _PRAD_GEM_TREE_H

#include <TClonesArray.h>
#include "GEMDataStruct.h"
#include "PRadDataHandler.h"
#include <vector>

#define CLONE_ARRAY
//#undef CLONE_ARRAY

class TTree;

class PRadGEMTree
{
public: 
    PRadGEMTree();
    ~PRadGEMTree();

    void Fill(const unsigned long &index, std::vector<GEMClusterStruct> &gem1, 
                                          std::vector<GEMClusterStruct> &gem2, 
                                          std::vector<HyCalHit> &hycal);
 
    void Save();
    void Clear();

    // public states
    float angular_resolution_from_moller1;
    float angular_resolution_from_moller2;
    float angle_moller1;
    float angle_moller2;
    float symm_dx;
    float symm_dy;
    float symm_coplanarity;

private:
    TFile *file;
    TTree *tree;

    unsigned long evt_id;
#ifdef CLONE_ARRAY
    TClonesArray *a_gem1;
    TClonesArray *a_gem2;
    TClonesArray *a_hycal;
#else   
    // different fill method
    int gem1_nhit;
    float gem1_x[100];
    float gem1_y[100];
    float gem1_x_charge[100];
    float gem1_y_charge[100];
    float gem1_energy[100];
    float gem1_z[100];
    int gem1_x_size[100];
    int gem1_y_size[100];
    
    int gem2_nhit;
    float gem2_x[100];
    float gem2_y[100];
    float gem2_x_charge[100];
    float gem2_y_charge[100];
    float gem2_energy[100];
    float gem2_z[100];
    int gem2_x_size[100];
    int gem2_y_size[100];
    
    int hycal_nhit;
    float hycal_x[100];
    float hycal_y[100];
    float hycal_energy[100];
#endif
};

#endif
