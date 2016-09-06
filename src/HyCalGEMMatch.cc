#include "HyCalGEMMatch.h"
#include "GEMCoord.h"
#include <TMath.h>
#include <iostream>

using namespace std;

HyCalGEMMatch::HyCalGEMMatch()
{
    delta = 60.;

    // mm +[ chamber + screw cap thickness]
    z_gem1 = 5300.;
    z_gem2 = 5260.;
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
    res_hycal_gem1.clear();
    res_hycal_gem2.clear();
}

void HyCalGEMMatch::Reset()
{
    Clear();
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

vector<HyCalHit> & HyCalGEMMatch::GetMatchHyCalGEM1()
{
    return res_hycal_gem1;
}

vector<HyCalHit> & HyCalGEMMatch::GetMatchHyCalGEM2()
{
    return res_hycal_gem2;
}

void HyCalGEMMatch::SetMatchCriteria(double & c)
{
    delta = c;
}

void HyCalGEMMatch::SetMatchCriteria(double && c)
{
    delta = c;
}

void HyCalGEMMatch::SetGEM1Z(double &z)
{
    z_gem1 = z;
}

void HyCalGEMMatch::SetGEM2Z(double &z)
{
    z_gem2 = z;
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
    else {
	//ShowMatch();
	return MetaMatch();
    }
}

int HyCalGEMMatch::MetaMatch()
{
    // for each hit on hycal, look for its matching point on gem
    double res = delta;
    int m_index;
    float m_e = 0.;
    bool match_gem1 = false;
    bool match_gem2 = false;

    int n_hycal = hycal_hit->size();
    for(int i=0;i<n_hycal;i++) 
    {
	match_gem1 = false;
	match_gem2 = false;
	res = delta;

	for(int ii=0;ii<2;ii++)
	{
	    int n_points = SinglePointMatchAlgorithm(res, gem[ii], hycal_hit, i, m_index);
	    if(n_points != 0)
	    {
		m_e = hycal_hit->at(i).E;
		if(ii == 0){
		    match_gem1 = true;
		    match_gem2 = false;
		}
		else {
		    match_gem1 = false;
		    match_gem2 = true;
		}
	    }
	}

	if(match_gem1 && match_gem2)
	    cout<<" HyCal GEM Match error..."
		<<endl;
	else if( match_gem1 ) {
	    gem[0][m_index].energy = m_e;
	    res_gem.push_back(gem[0][m_index]);
	    res_hycal.push_back(hycal_hit->at(i));
	}
	else if( match_gem2 ) {
	    gem[1][m_index].energy = m_e;
	    res_gem.push_back(gem[1][m_index]);
	    res_hycal.push_back(hycal_hit->at(i));
	}
    }
    return res_gem.size();
}

void HyCalGEMMatch::ProjectHyCalToPlaneZ(double &x, double &y, double &z)
{
    x = x*z/z_hycal;
    y = y*z/z_hycal;
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
    //SetMatchCriteria(10.); 
    // Hard code a stricter match criteria
    // b/c overlapping area is 44mm
    gem_coord -> GetPlaneClusterPlusMode(gem[0], gem[1]);
    if(gem[0].size() != 0)
	MetaMatchByGEM(gem[0], res_gem1, res_hycal_gem1);
    if(gem[1].size() != 0)
	MetaMatchByGEM(gem[1], res_gem2, res_hycal_gem2);
}

void HyCalGEMMatch::MetaMatchByGEM(vector<GEMClusterStruct> & _gem, vector<GEMClusterStruct> & _res, vector<HyCalHit> &_hycal)
{
    double res = delta;
    int m_index;
    float m_e = 0.;
    bool match_gem = false;

    int n_hycal = hycal_hit->size();
    int n_gem = _gem.size();

    double z = _gem[0].z;
    if( z<5200. || z > 5400)
	cout<<" HyCal GEM Match: Z ERROR."<<endl;

    for(int i=0;i<n_hycal;i++) 
    {
	res = delta;
	match_gem = false;
	int n_points = SinglePointMatchAlgorithm(res, _gem, hycal_hit, i, m_index);
	if( n_points != 0){
	    match_gem = true;
	    m_e = hycal_hit->at(i).E;
	}

	if( match_gem ) {
	    _gem[m_index].energy = m_e;
	    _res.push_back(_gem[m_index]);
	    _hycal.push_back(hycal_hit->at(i));
	}
    }
}

int HyCalGEMMatch::SinglePointMatchAlgorithm(double &match_criteria, vector<GEMClusterStruct> &gem, 
	vector<HyCalHit> * hycal, int & index, int &match_index)
{
    // requirement: gem only contain points from one chamber, ( same  z_gem )
    int points_within_range = 0;

    int n_gem = gem.size();
    double z = gem[0].z;
    if(z<5200. || z> 5400)
	cout<<" hycal gem match, z error..."<<endl;

    double x = (hycal_hit->at(index).x);
    double y = (hycal_hit->at(index).y);
    ProjectHyCalToPlaneZ(x, y, z);

    double dr = 0.;
    for(int j=0;j<n_gem;j++) 
    {
	dr = r(x-gem[j].x, y-gem[j].y);
	if(dr < match_criteria) 
	{
	    match_criteria = dr;
	    match_index = j;
	    points_within_range++;
	}
    }
    return points_within_range;
}

#include <TGraph.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TEllipse.h>

void HyCalGEMMatch::ShowMatch()
{
    if(gem[0].size() == 0 || gem[1].size() == 0)
	return;
    if(hycal_hit->size() == 0)
	return;

    int ng1 = gem[0].size();
    int ng2 = gem[1].size();
    int nh = hycal_hit->size();

    Float_t xg1[ng1], yg1[ng1];
    Float_t xg2[ng2], yg2[ng2];
    Float_t xh[nh], yh[nh];

    int index = 0;
    for(auto &i: gem[0]){
	xg1[index] = i.x;
	yg1[index++] = i.y;
    }
    index = 0;
    for(auto &i: gem[1]){
	xg2[index] = i.x;
	yg2[index++] = i.y;
    }
    index = 0;
    for(auto &i: *hycal_hit){
	xh[index] = i.x;
	xh[index++] = i.y;
    }

    TFile *f = new TFile("temp.root", "recreate");
    TCanvas *c = new TCanvas("c", "c", 800,800);
    TGraph *g1 = new TGraph(ng1, xg1, yg1);
    TGraph *g2 = new TGraph(ng2, xg2, yg2);
    TGraph *gh = new TGraph(nh, xh, yh);

    TEllipse *el[nh];
    for(int i=0;i<nh;i++){
	el[i] = new TEllipse(xh[i], yh[i], 60);
    }

    g1->SetName("g1");
    g2->SetName("g2");
    gh->SetName("gh");

    g1->SetMarkerStyle(22);
    g2->SetMarkerStyle(23);
    gh->SetMarkerStyle(8);

    g1->SetMarkerColor(2);
    g2->SetMarkerColor(4);
    gh->SetMarkerColor(6);

    g1->SetMarkerSize(1);
    g2->SetMarkerSize(1);
    gh->SetMarkerSize(1);
    c->DrawFrame(-600,-600,600,600);
    cout<<ng1<<", "<<ng2<<", "<<nh<<endl;
    for(int i=0;i<nh;i++){
	el[i]->SetLineColor(4);
	el[i]->Draw("P");
    }
    g1->Draw("P");
    g2->Draw("P");
    gh->Draw("P");
    c->Update();
    cout<<"press enter to continue..."<<endl;
    c->Write();
    f->Write();
    f->Close();
    getchar();
}
