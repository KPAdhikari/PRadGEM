//==============================================================================//
// A nasty class, dealing with Coordinate Transformation                        //
// Version 1:                                                                   //
//     All gem clusters will be projected to z=z_gem2 plane.                    //
// Version 2:                                                                   //
//     No project any more, All gem clusters in their own plane.                //
//     Version 1 is causing confusion.  In fact "z_gem1" and "z_gem2" are unknown.
//     Offsets also must be in their own Frame.                                 //
//     If need to convert them into one plane, make sure you have right offsts  //
//                                                                              //
// Xinzhan Bai                                                                  //
// 08/09/2016                                                                   //
//==============================================================================//
#include "GEMCoord.h"
#include "GEMOnlineHitDecoder.h"
#include "GEMCluster.h"
#include "GEMMapping.h"
#include <TH1F.h>

//#define PROJECT_TO_ONE_PLANE

using namespace std;

GEMCoord::GEMCoord()
{
    cout<<"GEMCoord Constructor..."
	<<endl;
    mapping = GEMMapping::GetInstance();

    InitPRadGeometry();
}

GEMCoord::~GEMCoord()
{
}

void GEMCoord::SetHitDecoder(GEMOnlineHitDecoder * p)
{
    hit_decoder = p;
}

void GEMCoord::GetClusterGEM( int nth, vector<GEMClusterStruct> &gem)
{
    // get charactoristics for nth gem detector
    //planeList.clear();
    fListOfClustersCleanFromPlane = hit_decoder->GetListOfClustersCleanFromPlane();
    if(fListOfClustersCleanFromPlane -> size() == 0)
        return;
    gem.clear();

    detector = mapping->GetDetectorFromID(nth);
    planeList = mapping->GetPlaneListFromDetector(detector);

    list<GEMCluster*> clusterX = (*fListOfClustersCleanFromPlane)[planeList.front()];
    list<GEMCluster*> clusterY = (*fListOfClustersCleanFromPlane)[planeList.back()];

    int sx = clusterX.size();
    int sy = clusterY.size();
    int N = (sx<sy)?sx:sy;

    list<GEMCluster*>::iterator itx = clusterX.begin();
    list<GEMCluster*>::iterator ity = clusterY.begin();

    for(int i=0;i<N;i++)
    {
	float c_x = (*itx)->GetClusterADCs();
	float c_y = (*ity)->GetClusterADCs();
	float x = (*itx)->GetClusterPosition();
	float y = (*ity)->GetClusterPosition();
	GEMClusterStruct cluster(x, y, c_x, c_y);
	cluster.x_size = (*itx)->GetNbOfHits();
	cluster.y_size = (*ity)->GetNbOfHits();
	cluster.chamber_id = nth;
	gem.push_back(cluster);
	itx++;
	ity++;
    }
}

void GEMCoord::GetClusterGEMPlusMode( int nth, vector<GEMClusterStruct> &gem)
{
    // get charactoristics for nth gem detector
    //planeList.clear();
    fListOfClustersCleanFromPlane = hit_decoder->GetListOfClustersCleanFromPlane();
    if(fListOfClustersCleanFromPlane -> size() == 0)
        return;
    gem.clear();

    detector = mapping->GetDetectorFromID(nth);
    planeList = mapping->GetPlaneListFromDetector(detector);

    list<GEMCluster*> clusterX = (*fListOfClustersCleanFromPlane)[planeList.front()];
    list<GEMCluster*> clusterY = (*fListOfClustersCleanFromPlane)[planeList.back()];

    list<GEMCluster*>::iterator itx = clusterX.begin();
    list<GEMCluster*>::iterator ity = clusterY.begin();

    for(itx=clusterX.begin();itx!=clusterX.end();++itx)
    {
	for(ity=clusterY.begin();ity!=clusterY.end();++ity)
	{
	    float c_x = (*itx)->GetClusterADCs();
	    float c_y = (*ity)->GetClusterADCs();
	    float x = (*itx)->GetClusterPosition();
	    float y = (*ity)->GetClusterPosition();
	    GEMClusterStruct cluster(x, y, c_x, c_y);
	    cluster.x_size = (*itx)->GetNbOfHits();
	    cluster.y_size = (*ity)->GetNbOfHits();
	    cluster.chamber_id = nth;
	    gem.push_back(cluster);
	}
    }
}

void GEMCoord::InitPRadGeometry()
{
    //--------------------------------------------
    // Convert coordinate
    // -------------------------------------------
    //    Move GEM coordinate to GEM Hole center
    //    overlapping area: hole diameter (44mm)
    //    origin shift: 550.4/2 - (44-pitch)/2 = 253.2
    //
    //    gem1: x = x-253.2; y = y
    //    gem2: x = -x+253.2; y=-y
    //    
    //    gem1 is the one away from target, 
    //    left one, if looking downstream
    //    
    //    right-hand coordinate, HyCal Y axis
    //    must be pointing downward
    //
    //    beam downstream is z axis direction
    //    chamber thickness + screw thickness = 20mm
    //    5304 and 5264 do not include this 20mm
    //---------------------------------------------

    origin_transfer = 253.2;
    overlap_length = 44.;
    z_gem[0] = 5304.;
    z_gem[1] = 5264.;

    gem_offset_x = 0.;
    gem_offset_y = 0.;
    gem_hycal_offset_x = 0.;
    gem_hycal_offset_y = 0.;
    gem_beam_offset_x = 0.;
    gem_beam_offset_y = 0.;
}

void GEMCoord::SetGEMOffsetX(double v)
{
    // gem_offset_x will be used in 
    // gem local frame

    // make sure the offset provided
    // is in gem1 local frame.

    gem_offset_x = v;
}

void GEMCoord::SetGEMOffsetY(double v)
{
    // gem_offset_y will be used in 
    // gem local frame

    // make sure the offset provided
    // is in gem1 local frame.

    gem_offset_y = v;
}

void GEMCoord::SetGEMHyCalOffsetX(double v)
{
    gem_hycal_offset_x = v;
}

void GEMCoord::SetGEMHyCalOffsetY(double v)
{
    gem_hycal_offset_y = v;
}

void GEMCoord::SetGEMBeamOffsetX(double v)
{
    gem_beam_offset_x = v;
}

void GEMCoord::SetGEMBeamOffsetY(double v)
{
    gem_beam_offset_y = v;
}

pair<double, double> GEMCoord::TransformGEMCoord(int n, pair<double, double> && c)
{
    //----------------------------------------
    //  Process:
    //      step 1), align gem1 to gem2, make 
    //               a gem plane based on gem2
    //      step 2), align gem plane to hycal
    //      step 3), align hycal to beam 
    // 
    //  Better: 
    //      hycal has worse position resolution
    //      align gem plane directly to beam
    //      put gem_hycal_offsets to zero
    //      so:
    //      step 1), align gem1 to gem2
    //      step 2), align gem plane to beam
    //----------------------------------------

    pair<double, double> res;

    //----------------------------------------
    //Meaning:
    //
    //  gem_offset_x, y: gem1 origin coord in
    //      gem2 local frame
    //
    //  gem_beam_offset_x, y: beam position coord
    //      in gem2 local frame
    //
    //  gem_hycal_offset_x, y: not useful.
    //----------------------------------------
    if(n==0)
    {
	res.first = (c.first - origin_transfer) + gem_offset_x - gem_hycal_offset_x - gem_beam_offset_x;
	res.second = -c.second + gem_offset_y - gem_hycal_offset_y - gem_beam_offset_y;
    }
    else if(n == 1)
    {
	res.first = (origin_transfer - c.first) - gem_hycal_offset_x - gem_beam_offset_x;
	res.second = c.second - gem_hycal_offset_y - gem_beam_offset_y;
    }
    else
	cout<<"Transform GEM Coord Error..."
	    <<endl;
    return res;
}

void GEMCoord::ProjectToGEM2(pair<double, double> & res)
{
    res.first = res.first * z_gem[1]/z_gem[0];
    res.second = res.second * z_gem[1]/z_gem[0];
}

void GEMCoord::GetPlaneClusterTransferred(vector<GEMClusterStruct> &gem1,
	vector<GEMClusterStruct> &gem2)
{
    /* keep all overlapping area*/
    fListOfClustersCleanFromPlane = hit_decoder->GetListOfClustersCleanFromPlane();
    list<GEMCluster*> cluster_x1 = (*fListOfClustersCleanFromPlane)["pRadGEM1X"];
    list<GEMCluster*> cluster_y1 = (*fListOfClustersCleanFromPlane)["pRadGEM1Y"];
    list<GEMCluster*> cluster_x2 = (*fListOfClustersCleanFromPlane)["pRadGEM2X"];
    list<GEMCluster*> cluster_y2 = (*fListOfClustersCleanFromPlane)["pRadGEM2Y"];

    XYClusterMatchNormalMode(0, gem1, cluster_x1, cluster_y1);
    XYClusterMatchNormalMode(1, gem2, cluster_x2, cluster_y2);
}

void GEMCoord::GetPlaneClusterCutMode(vector<GEMClusterStruct> &gem1,
	vector<GEMClusterStruct> &gem2)
{
    /* cut overlapping area */
    fListOfClustersCleanFromPlane = hit_decoder->GetListOfClustersCleanFromPlane();
    list<GEMCluster*> cluster_x1 = (*fListOfClustersCleanFromPlane)["pRadGEM1X"];
    list<GEMCluster*> cluster_y1 = (*fListOfClustersCleanFromPlane)["pRadGEM1Y"];
    list<GEMCluster*> cluster_x2 = (*fListOfClustersCleanFromPlane)["pRadGEM2X"];
    list<GEMCluster*> cluster_y2 = (*fListOfClustersCleanFromPlane)["pRadGEM2Y"];

    XYClusterMatchCutMode(0, gem1, cluster_x1, cluster_y1);
    XYClusterMatchCutMode(1, gem2, cluster_x2, cluster_y2);
}

void GEMCoord::GetPlaneClusterPlusMode(vector<GEMClusterStruct> &gem1,
	vector<GEMClusterStruct> &gem2)
{
    fListOfClustersCleanFromPlane = hit_decoder->GetListOfClustersCleanFromPlane();

    list<GEMCluster*> cluster_x1 = (*fListOfClustersCleanFromPlane)["pRadGEM1X"];
    list<GEMCluster*> cluster_y1 = (*fListOfClustersCleanFromPlane)["pRadGEM1Y"];
    list<GEMCluster*> cluster_x2 = (*fListOfClustersCleanFromPlane)["pRadGEM2X"];
    list<GEMCluster*> cluster_y2 = (*fListOfClustersCleanFromPlane)["pRadGEM2Y"];

    XYClusterMatchPlusMode(0, gem1, cluster_x1, cluster_y1);
    XYClusterMatchPlusMode(1, gem2, cluster_x2, cluster_y2);
}

void GEMCoord::XYClusterMatchPlusMode(int detID, vector<GEMClusterStruct> &gem, 
	list<GEMCluster*> &xlist, list<GEMCluster*> &ylist)
{
    for( auto &i : xlist)
    {
	for(auto &j : ylist)
	{
	    pair<double, double> res;
	    res = TransformGEMCoord(detID, make_pair(i->GetClusterPosition(), j->GetClusterPosition()));
#ifdef PROJECT_TO_ONE_PLANE
	    if(detID == 0)
		ProjectToGEM2(res); // project to gem2
#endif
	    float c_x = i->GetClusterADCs();
	    float c_y = j->GetClusterADCs();
	    float x = res.first;
	    float y = res.second;
	    GEMClusterStruct cluster(x,y,c_x,c_y);
	    cluster.x_size = i->GetNbOfHits();
	    cluster.y_size = j->GetNbOfHits();
	    cluster.chamber_id = detID;
	    cluster.z = z_gem[detID];
	    gem.push_back(cluster);
	}
    }
}

void GEMCoord::XYClusterMatchCutMode(int detID, vector<GEMClusterStruct> &gem, 
	list<GEMCluster*> &xlist, list<GEMCluster*> &ylist)
{
    int s1 = xlist.size();
    int s2 = ylist.size();  
    int nbCluster = (s1<s2)?s1:s2;

    // cutting edge
    float edge = 0;

    if(nbCluster != 0)
    {
	list<GEMCluster*>::iterator itx = xlist.begin();
	list<GEMCluster*>::iterator ity = ylist.begin();
	for(int i = 0;i<nbCluster;i++)
	{
	    pair<double, double> res;
	    res = TransformGEMCoord(detID, make_pair( (*itx)->GetClusterPosition(),(*ity)->GetClusterPosition()));
#ifdef PROJECT_TO_ONE_PLANE
	    if( detID == 0)
		ProjectToGEM2(res); // project to gem2
#endif
	    if( (detID == 0 && res.first <= edge) || (detID == 1 && res.first >= edge)) 
	    {
		float c_x = (*itx)->GetClusterADCs();
		float c_y = (*ity)->GetClusterADCs();
		float x = res.first;
		float y = res.second;
		GEMClusterStruct cluster(x, y, c_x, c_y);
		cluster.x_size = (*itx)->GetNbOfHits();
		cluster.y_size = (*ity)->GetNbOfHits();
		cluster.chamber_id = detID;
		cluster.z = z_gem[detID];
		gem.push_back(cluster);
		*itx++;
		*ity++;
	    }
	}
    }
}

void GEMCoord::XYClusterMatchNormalMode(int detID, vector<GEMClusterStruct> &gem, 
	list<GEMCluster*> &xlist, list<GEMCluster*> &ylist)
{
    int s1 = xlist.size();
    int s2 = ylist.size();  
    int nbCluster = (s1<s2)?s1:s2;

    if(nbCluster != 0)
    {
	list<GEMCluster*>::iterator itx = xlist.begin();
	list<GEMCluster*>::iterator ity = ylist.begin();
	for(int i = 0;i<nbCluster;i++)
	{
	    pair<double, double> res;
	    res = TransformGEMCoord(detID, make_pair((*itx)->GetClusterPosition(),(*ity)->GetClusterPosition()));
#ifdef PROJECT_TO_ONE_PLANE
	    if( detID == 0)
		ProjectToGEM2(res); // project to gem2
#endif
	    float c_x = (*itx)->GetClusterADCs();
	    float c_y = (*ity)->GetClusterADCs();
	    float x = res.first;
	    float y = res.second;
	    GEMClusterStruct cluster(x, y, c_x, c_y);
	    cluster.x_size = (*itx)->GetNbOfHits();
	    cluster.y_size = (*ity)->GetNbOfHits();
	    cluster.chamber_id = detID;
	    cluster.z = z_gem[detID];
	    gem.push_back(cluster);
	    *itx++;
	    *ity++;
	}
    }
}


void GEMCoord::FillHistos(TH1F* hNbClusterPerPlaneX[], 
	TH1F* hNbClusterPerPlaneY[], 
	TH1F* hClusterDistX[], 
	TH1F* hClusterDistY[])
{
    fListOfClustersCleanFromPlane = hit_decoder->GetListOfClustersCleanFromPlane();

    list<GEMCluster*> cluster_x1 = (*fListOfClustersCleanFromPlane)["pRadGEM1X"];
    list<GEMCluster*> cluster_y1 = (*fListOfClustersCleanFromPlane)["pRadGEM1Y"];
    list<GEMCluster*> cluster_x2 = (*fListOfClustersCleanFromPlane)["pRadGEM2X"];
    list<GEMCluster*> cluster_y2 = (*fListOfClustersCleanFromPlane)["pRadGEM2Y"];

    int s1 = cluster_x1.size();
    hNbClusterPerPlaneX[0]->Fill(s1);
    int s2 = cluster_y1.size();
    hNbClusterPerPlaneY[0]->Fill(s2);
    for(int i=0;i<s1;i++)
    {
	list<GEMCluster*>::iterator itx = cluster_x1.begin();
	hClusterDistX[0] -> Fill( ( *itx++)->GetClusterPosition() );
    }
    for(int i=0;i<s2;i++)
    {
	list<GEMCluster*>::iterator ity = cluster_y1.begin();
	hClusterDistY[0]->Fill( (*ity++)->GetClusterPosition() );
    }

    s1 = cluster_x2.size();
    hNbClusterPerPlaneX[1] -> Fill(s1);
    s2 = cluster_y2.size();
    hNbClusterPerPlaneY[1] -> Fill(s2);
    for(int i=0;i<s1;i++)
    {
	list<GEMCluster*>::iterator itx = cluster_x2.begin();
	hClusterDistX[1] -> Fill( ( *itx++)->GetClusterPosition() );
    }
    for(int i=0;i<s2;i++)
    {
	list<GEMCluster*>::iterator ity = cluster_y2.begin();
	hClusterDistY[1]->Fill( (*ity++)->GetClusterPosition() );
    }
}
