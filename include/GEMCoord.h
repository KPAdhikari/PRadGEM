#ifndef GEMCOORD_H
#define GEMCOORD_H
#include <list>
#include "GEMDataStruct.h"
#include <vector>
#include <map>
#include <TString.h>

class GEMCluster;
class GEMMapping;
class GEMOnlineHitDecoder;
class TH1F;

class GEMCoord
{
public:
    GEMCoord();
    ~GEMCoord();
    void SetHitDecoder(GEMOnlineHitDecoder *);

    //get charactoristics for nth gem detector
    void GetClusterGEM(int & nth, std::vector<GEMClusterStruct> &gem);

private:
    GEMOnlineHitDecoder * hit_decoder;
    GEMMapping * mapping;
    std::map<TString, std::list<GEMCluster*> > * fListOfClustersCleanFromPlane;

    int nDet;
    TString detector;
    std::list<TString> planeList;

    /* following only for PRad */
public:
    void GetPlaneClusterTransferred(std::vector<GEMClusterStruct> &gem1,
	    std::vector<GEMClusterStruct> &gem2);
    void GetPlaneClusterCutMode(std::vector<GEMClusterStruct> &gem1,
	    std::vector<GEMClusterStruct> &gem2);
    void GetPlaneClusterPlusMode(std::vector<GEMClusterStruct> &gem1,
	    std::vector<GEMClusterStruct> &gem2);

    void XYClusterMatchNormalMode(int , std::vector<GEMClusterStruct> &, 
	    std::list<GEMCluster*> &, std::list<GEMCluster*> &);
    void XYClusterMatchCutMode(int , std::vector<GEMClusterStruct> &, 
	    std::list<GEMCluster*> &, std::list<GEMCluster*> &);
    void XYClusterMatchPlusMode(int , std::vector<GEMClusterStruct> &, 
	    std::list<GEMCluster*> &, std::list<GEMCluster*> &);

    void FillHistos(TH1F**, TH1F**, TH1F**, TH1F**);

    void InitPRadGeometry();
    std::pair<double, double> TransformGEMCoord(int, std::pair<double, double> &&);
    void ProjectToGEM2(std::pair<double, double> &);
    void SetGEMOffsetX(double );
    void SetGEMOffsetY(double );
    void SetGEMHyCalOffsetX(double );
    void SetGEMHyCalOffsetY(double );
    void SetGEMBeamOffsetX(double );
    void SetGEMBeamOffsetY(double );

private:
    double origin_transfer;
    double overlap_length;
    double z_gem[2];

    double gem_offset_x;
    double gem_offset_y;
    double gem_hycal_offset_x;
    double gem_hycal_offset_y;
    double gem_beam_offset_x;
    double gem_beam_offset_y;
};

#endif
