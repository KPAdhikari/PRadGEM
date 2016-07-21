#ifndef _PRAD_GEM_TREE_H
#define _PRAD_GEM_TREE_H

#include <TClonesArray.h>
#include "GEMDataStruct.h"
#include "PRadDataHandler.h"
#include <vector>

class TTree;

class PRadGEMTree
{
public: 
    PRadGEMTree();
    ~PRadGEMTree();

    void Fill(const unsigned long &index, const std::vector<GEMClusterStruct> &gem1, 
                                          const std::vector<GEMClusterStruct> &gem2, 
                                          const std::vector<HyCalHit> &hycal);
 
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
    TTree *tree;

    unsigned long evt_id;

    TClonesArray *a_gem1;
    TClonesArray *a_gem2;
    TClonesArray *a_hycal;

};

#endif
