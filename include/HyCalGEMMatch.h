#ifndef HYCALGEMMATCH_H
#define HYCALGEMMATCH_H
#include <vector>
#include <cstring>
#include "datastruct.h"
#include "GEMDataStruct.h"
#include "PRadEventStruct.h"

class GEMCoord;

class HyCalGEMMatch
{
public:
    HyCalGEMMatch();
    ~HyCalGEMMatch();
    void SetGEMCoord(GEMCoord * fgemcoord);
    void SetHyCalVector(std::vector<HyCalHit> * hycalhit);
    std::vector<GEMClusterStruct> & GetMatchGEM();
    std::vector<GEMClusterStruct> & GetMatchGEM1();
    std::vector<GEMClusterStruct> & GetMatchGEM2();
    std::vector<HyCalHit> & GetMatchHyCal();
    void Clear();

    // Note 1:
    // match starting from hycal, for each cluster on hycal
    // look up matching cluster on gem. 
    // gem matching results will be organized into one plane.
    int Match();
    int MetaMatch();
    // Note 2: 
    // match starting from hycal, for each cluster on hycal
    // look up matching cluster on gem.
    // for single-gem moller events to detect offset
    // each gem chamber has its own matching result.
    void MatchByGEM(); 
    void MetaMatchByGEM(std::vector<GEMClusterStruct> & _gem, 
	    std::vector<GEMClusterStruct> & _res);
    // Note 3:
    // match starting from gem, for each cluster on gem
    // look up matching cluster on hycal.
    // for calibration run, distinguish e+ from e-
    // each gem gem chamber has its own matching result.
    //void MatchAccordToGEM();
    //void MetaMatchAccordToGEM(std::vector<GEMClusterStruct> & _gem, 
//	    std::vector<GEMClusterStruct> & _res);

    void SetMatchCriteria(double & );
    void SetMatchCriteria(double &&);
    void SetGEMZ(double &);
    void SetHyCalZ(double &);
    double r(double &, double &);
    double r(double &&, double &&);
    void ProjectHyCalToGEM(double &, double &);

private:
    GEMCoord *gem_coord;
    std::vector<GEMClusterStruct> res_gem;
    std::vector<GEMClusterStruct> res_gem1;
    std::vector<GEMClusterStruct> res_gem2;
    std::vector<HyCalHit> res_hycal;
    std::vector<GEMClusterStruct> gem[2];
    std::vector<HyCalHit> *hycal_hit;

    double z_gem;
    double z_hycal;
    double delta;
};

#endif
