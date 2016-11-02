#ifndef GEMTREE_H
#define GEMTREE_H
#include "GEMDataStruct.h"
#include <vector>
#include <unordered_map>
#include <cstring>
#include "datastruct.h"
#include "PRadEventStruct.h"
#include <unordered_map>

#define NDETECTOR 2

class TTree;
class TFile;
class TH1F;
class PRadMoller;
class PRadEP;
class GEMConfigure;

class GEMTree
{
    // ------------------------------- common setting --------------------------------------
public:
    GEMTree();
    ~GEMTree();
    void WriteToDisk();
    void SetGEMConfigure(GEMConfigure *con);

private:
    TFile *file;
    GEMConfigure *configure;

    // for all trees
    unsigned int evt_id;


    // ------------------------------- gem tree -------------------------------------
public:
    void PushCalibrationData(std::vector<GEMClusterStruct> & gem1,
	    std::vector<GEMClusterStruct> &gem2,
	    std::unordered_multimap<std::string, double> &tdc_map,
	    std::vector<HyCalHit> *hycalhit);

    void PushDetector(int, std::vector<GEMClusterStruct>);
    void PushTDCValue(std::unordered_multimap<std::string, double> &tdc_map);
    void PushHyCalData(std::vector<HyCalHit> *hycalhit);

    void FillGEMTree();
    void InitGEMTree(int ndet);
    void InitTDCGroup();

private:
    TTree *gem_tree;

    bool empty_gem_event;
    // number of hits on each detector
    const int static _nhits = 2000; // if _nhits is not large enough, it will become a potential break threat
    std::string scin_group[2];
    int hycal_group_q;
    std::string hycal_group[9];

    // max 4 detectors
    int nhits[NDETECTOR];
    float x[NDETECTOR][_nhits];
    float y[NDETECTOR][_nhits];
    float x_charge[NDETECTOR][_nhits];
    float y_charge[NDETECTOR][_nhits];
    float energy[NDETECTOR][_nhits];
    float z[NDETECTOR][_nhits];
    float x_size[NDETECTOR][_nhits];
    float y_size[NDETECTOR][_nhits];

    // tdc information
    int n_S2;
    int n_S1;
    int n_hycal;
    float TDCS2[_nhits];
    float TDCS1[_nhits];
    float TDCHyCal[_nhits];

    // hycal information
    int hycal_nhits;
    float hycal_x[_nhits];
    float hycal_y[_nhits];
    float hycal_e[_nhits];


    // ------------------------------- epics tree -------------------------------------
public:
    void PushEpics(std::unordered_map<std::string, double> &);
    void FillEpicsTree();
    void InitEpicsTree();
    void SetEventID(unsigned int &id);

private:
    TTree *epics_tree;
    bool empty_epics_event;
    EpicsStruct _epics_tree;
    int event_id;


    // ------------------------------- physics tree -------------------------------------
public:
    void InitPhysicsTree();
    void PushMoller(PRadMoller *);
    void PushEP(PRadEP *);
    void FillMollerTree();
    void FillEPTree();

private:
    // moller
    bool empty_moller_event;
    TTree *moller_tree;
    moller_data_t moller_data;
    float moller_center_x;
    float moller_center_y;
    float moller_pos_res_dx;
    float moller_pos_res_dy;
    // ep 
    bool empty_ep_event;
    TTree *ep_tree;
    ep_data_t ep_data;
    float q_square;
    // ep + moller
    TTree *ep_moller_tree;
    const int static _ncluster = 10;
    int nCluster;
    float scatt_energy[_ncluster];
    float scatt_angle[_ncluster];
    float scatt_x[_ncluster];
    float scatt_y[_ncluster];
    int detector_id[_ncluster];


    // ----------------------- calibration data offset tree -----------------------------
public:
    void InitCaliOffsetTree();
    void PushCaliOffset(double, double);
    void FillCaliOffsetTree();

private:
    TTree *cali_offset_tree;
    double cali_x_offset;
    double cali_y_offset;


    // ------------------------- production data offset tree -----------------------------
public:
    void InitProdOffsetTree();
    void PushProdOffset(int , PRadMoller *);
    void FillProdOffsetTree();

private:
    // offset from production runs, overlapping area method
    // gem1
    TTree *prod_offset_tree1;
    double prod_gem1_dx;
    double prod_gem1_dy;
    double prod_gem1_energy1;
    double prod_gem1_energy2;
    double prod_gem1_angle1;
    double prod_gem1_angle2;
    double prod_gem1_coplanarity;
    double prod_gem1_open_angle;
    int prod_gem1_ncluster;
    double prod_gem1_scattx[_ncluster];
    double prod_gem1_scatty[_ncluster];
    bool empty_prod_gem1;
    // gem2
    TTree *prod_offset_tree2;
    double prod_gem2_dx;
    double prod_gem2_dy;
    double prod_gem2_energy1;
    double prod_gem2_energy2;
    double prod_gem2_angle1;
    double prod_gem2_angle2;
    double prod_gem2_coplanarity;
    double prod_gem2_open_angle;
    int prod_gem2_ncluster;
    double prod_gem2_scattx[_ncluster];
    double prod_gem2_scatty[_ncluster];
    bool empty_prod_gem2;
    // overall offset
    TTree *prod_offset_tree_res;
    double prod_offset_x;
    double prod_offset_y;


    // --------------------------- overlapp area data tree -------------------------------
public:
    void InitOverlapTree();
    void PushOverlapMollerTree(PRadMoller*, PRadMoller*);
    void PushOverlapEpTree(PRadEP*, PRadEP*);
    void FillOverlapTree();

private:
    bool olp_empty_moller_event;
    TTree *overlap_moller_tree;
    bool olp_empty_ep_event;
    TTree *overlap_ep_tree;

    moller_data_t olp_moller_data1;
    moller_data_t olp_moller_data2;
    ep_data_t olp_ep_data1;
    ep_data_t olp_ep_data2;
};

#endif
