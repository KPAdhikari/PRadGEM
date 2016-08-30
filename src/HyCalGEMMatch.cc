#include "HyCalGEMMatch.h"
#include "GEMCoord.h"
#include <TMath.h>
#include <iostream>

using namespace std;

HyCalGEMMatch::HyCalGEMMatch()
{
    delta = 30.;

    // mm +[ chamber + screw cap thickness]
    // in gem_coord, coords already projected 
    // to gem2
    z_gem = 5260. + 20.;
    z_hycal = 5600.;
}

HyCalGEMMatch::~HyCalGEMMatch()
{
    Clear();
}

void HyCalGEMMatch::SetGEMCoord(GEMCoord *fgemcoord)
{
    gem_coord = fgemcoord;
}

void HyCalGEMMatch::SetHyCalVector(vector<HyCalHit> *hycalhit)
{
    hycal_hit = hycalhit;
}

void HyCalGEMMatch::Clear()
{
    gem[0].clear();
    gem[1].clear();
    res_gem.clear();
    res_gem1.clear();
    res_gem2.clear();
    res_hycal.clear();
}

vector<GEMClusterStruct> & HyCalGEMMatch::GetMatchGEM()
{
    return res_gem;
}

vector<GEMClusterStruct> & HyCalGEMMatch::GetMatchGEM1()
{
    return res_gem1;
}

vector<GEMClusterStruct> & HyCalGEMMatch::GetMatchGEM2()
{
    return res_gem2;
}

vector<HyCalHit> & HyCalGEMMatch::GetMatchHyCal()
{
    return res_hycal;
}

void HyCalGEMMatch::SetMatchCriteria(double & c)
{
    delta = c;
}

void HyCalGEMMatch::SetMatchCriteria(double && c)
{
    delta = c;
}

void HyCalGEMMatch::SetGEMZ(double &z)
{
    z_gem = z;
}

void HyCalGEMMatch::SetHyCalZ(double & z)
{
    z_hycal = z;
}

int HyCalGEMMatch::Match()
{
    Clear();
    gem_coord -> GetPlaneClusterPlusMode(gem[0], gem[1]);
    if( (gem[0].size() ==0) && (gem[1].size() == 0) )
        return 0;
    else if( hycal_hit -> size() == 0)
        return 0;
    else 
        return MetaMatch();
}

int HyCalGEMMatch::MetaMatch()
{
    // for each hit on hycal, look for its matching point on gem
    double res = delta;
    double dr = 0.;
    int m_index;
    float m_e = 0.;
    bool match_gem1 = false;
    bool match_gem2 = false;

    int n_hycal = hycal_hit->size();
    for(int i=0;i<n_hycal;i++) 
    {
	double x = (hycal_hit->at(i).x);
	double y = (hycal_hit->at(i).y);
	ProjectHyCalToGEM(x, y);
	match_gem1 = false;
	match_gem2 = false;
	res = delta;

	for(int ii=0;ii<2;ii++){
	    for(int j=0;j<gem[ii].size();j++){
		dr = r(x-gem[ii][j].x, y-gem[ii][j].y);
		if(dr < res) {
		    res = dr;
		    m_index = j;
		    m_e = hycal_hit->at(i).E;
		    if(ii == 0){
			match_gem1 = true;
			match_gem2 = false;
		    }
		    else{
			match_gem2 = true;
			match_gem1 = false;
		    }
		}
	    }
	}
	if(match_gem1 && match_gem2)
	    cout<<" HyCal GEM Match error..."
		<<endl;
	else if( match_gem1 ) {
	    gem[0][m_index].energy = m_e;
	    gem[0][m_index].z = z_gem;
	    res_gem.push_back(gem[0][m_index]);
	    res_hycal.push_back(hycal_hit->at(i));
	}
	else if( match_gem2 ) {
	    gem[1][m_index].energy = m_e;
	    gem[1][m_index].z = z_gem;
	    res_gem.push_back(gem[1][m_index]);
	    res_hycal.push_back(hycal_hit->at(i));
	}
    }
    return res_gem.size();
}

void HyCalGEMMatch::ProjectHyCalToGEM(double &x, double &y)
{
    x = x*z_gem/z_hycal;
    y = y*z_gem/z_hycal;
}

double HyCalGEMMatch::r(double & x, double & y)
{
    return TMath::Sqrt(x*x + y*y);
}

double HyCalGEMMatch::r(double && x, double && y)
{
    return TMath::Sqrt(x*x + y*y);
}

void HyCalGEMMatch::MatchByGEM()
{
    Clear();
    SetMatchCriteria(10.); 
    // Hard code a stricter match criteria
    // b/c overlapping area is 44mm
    gem_coord -> GetPlaneClusterPlusMode(gem[0], gem[1]);
    if(gem[0].size() != 0)
	MetaMatchByGEM(gem[0], res_gem1);
    if(gem[1].size() != 0)
	MetaMatchByGEM(gem[1], res_gem2);
}

void HyCalGEMMatch::MetaMatchByGEM(vector<GEMClusterStruct> & _gem, vector<GEMClusterStruct> & _res)
{
    double res = delta;
    double dr = 0.;
    int m_index;
    float m_e = 0.;
    bool match_gem = false;

    int n_hycal = hycal_hit->size();
    int n_gem = _gem.size();
    for(int i=0;i<n_hycal;i++) {
	double x = (hycal_hit->at(i).x);
	double y = (hycal_hit->at(i).y);
	ProjectHyCalToGEM(x, y);

	match_gem = false;
	res = delta;
	for(int j=0;j<n_gem;j++) {
	    dr = r(x-_gem[j].x, y-_gem[j].y);
	    if(dr < res) {
		res = dr;
		m_index = j;
		m_e = hycal_hit->at(i).E;
		match_gem = true;
	    }
	}
	if( match_gem ) {
	    _gem[m_index].energy = m_e;
	    _gem[m_index].z = z_gem;
	    _res.push_back(_gem[m_index]);
	}
    }
}
