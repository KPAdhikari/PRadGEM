#ifndef GEMTREE_H
#define GEMTREE_H
#include "GEMDataStruct.h"
#include <vector>
#include <cstring>
#include "datastruct.h"
#include "PRadEventStruct.h"
#include <unordered_map>

#define NDETECTOR 4

class TTree;
class TFile;
class TH1F;
class PRadMoller;
class PRadEP;

class GEMTree
{
public:
    struct moller_data_t
    {
        unsigned int event_id;
	int chamber_id1;
        float x1;
	float y1;
	float x1_hycal;
	float y1_hycal;
	float e1;
	float angle1;
	int chamber_id2;
	float x2;
	float y2;
	float x2_hycal;
	float y2_hycal;
	float e2;
	float angle2;

	float coplanarity;

	moller_data_t()
        : event_id(0), chamber_id1(0),  x1(0),y1(0), 
	x1_hycal(0), y1_hycal(0),
	e1(0), angle1(0), chamber_id2(0), x2(0), y2(0), 
	x2_hycal(0), y2_hycal(0),
	e2(0), angle2(0), coplanarity(0)
	{}
    };

    struct ep_data_t
    {
        unsigned int event_id;
	int chamber_id;
	float x;
	float y;
	float x_hycal;
	float y_hycal;
	float e;
	float angle;

	ep_data_t()
	: event_id(0), chamber_id(0), x(0),
	y(0), x_hycal(0), y_hycal(0), e(0), angle(0)
	{}
    };

    GEMTree();
    ~GEMTree();
    // GEM TREE
    void PushDetector(int, std::vector<GEMClusterStruct>);
    void FillGEMTree();
    void InitGEMTree(int ndet);

    void WriteToDisk();

    //EPICS TREE
    void PushEpics(std::unordered_map<std::string, double> &);
    void FillEpicsTree();
    void InitEpicsTree();
    void SetEventID(unsigned int &id);

    // PHYSICS TREE
    struct Energy_Angle
    {
        double energy;
	double angle;
	Energy_Angle(double e, double a)
	: energy(e), angle(a)
	{}
	Energy_Angle()
	: energy(-10000), angle(-10000)
	{}
    };
    void InitPhysicsTree();
    void PushMoller(PRadMoller *);
    void PushEP(PRadEP *);
    void FillMollerTree();
    void FillEPTree();

    // calibration offset tree
    void InitCaliOffsetTree();
    void PushCaliOffset(double, double);
    void FillCaliOffsetTree();

    // production offset tree
    void InitProdOffsetTree();
    void PushProdOffset(int , PRadMoller *);
    void FillProdOffsetTree();

    // overlap area tree
    // ep tree
    void InitOverlapTree();
    void PushOverlapMollerTree(PRadMoller*, PRadMoller*);
    void PushOverlapEpTree(PRadEP*, PRadEP*);
    void FillOverlapTree();

    // moller tree

    // can also save histograms
    // make histos public for easy access
    TH1F *adc_x1;

private:
    TFile *file;

    // for all trees
    unsigned int evt_id;

    // GEM TREE
    TTree *gem_tree;
    int nDet;
    bool empty_gem_event;
    // number of hits on each detector
    const int static _nhits = 100;
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

    // EPICS TREE
    TTree *epics_tree;
    bool empty_epics_event;
    EpicsStruct _epics_tree;
    int event_id;

    //PHYSICS TREE
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

    // offset from calibration runs
    TTree *cali_offset_tree;
    double cali_x_offset;
    double cali_y_offset;

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

    // overlap area tree
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
