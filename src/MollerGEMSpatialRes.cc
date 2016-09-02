#include "MollerGEMSpatialRes.h"
#include <TMath.h>

#define RES_UNDEFINED 9999
#define PI 3.1415926

using namespace std;

MollerGEMSpatialRes::MollerGEMSpatialRes()
    : beam_energy(2147.), x_origin(0.),
      y_origin(0.), z_gem1(5300.+20.), z_gem2(5260.+20.)
{
    diff.first = RES_UNDEFINED;
    diff.second = RES_UNDEFINED;
}

MollerGEMSpatialRes::~MollerGEMSpatialRes()
{
}

void MollerGEMSpatialRes::Process(vector<pair<int, pair<double, double> > >in)
{
    Reset();

    // from "far" map "near" leads to a smaller diff
    // from "near" map "far" leads to a larger diff
    // which one to use?
    if( r(in[1].second.first, in[1].second.second) < r(in[0].second.first, in[0].second.second) ){
	SetPoint1(in[0].second);
	SetPoint2(in[1].second);
	if( in[0].first == 0){
	    z_p1 = z_gem1;
	    z_p2 = z_gem2;
	}
	else{
	    z_p1 = z_gem2;
	    z_p2 = z_gem1;
	}
    }
    else{
	SetPoint1(in[1].second);
	SetPoint2(in[0].second);
	if( in[1].first == 0){
	    z_p1 = z_gem1;
	    z_p2 = z_gem2;
	}
	else{
	    z_p1 = z_gem2;
	    z_p2 = z_gem1;
	}
    }

    p1.first -= x_origin;
    p1.second -= y_origin;
    p2.first -= x_origin;
    p2.second -= y_origin;
    _Process();
}

void MollerGEMSpatialRes::Reset()
{
    diff.first = RES_UNDEFINED;
    diff.second = RES_UNDEFINED;
}

void MollerGEMSpatialRes::_Process()
{
    double r_1 = r(p1.first, p1.second);
    double e_1 = EnergyFromAngle(TMath::ATan(r_1/z_p1));
    double e_2 = beam_energy - e_1; 
    double r_2 = TMath::Tan( AngleFromEnergy(e_2) )* z_p2;

    double x_t = - p1.first * r_2/r_1;
    double y_t = -p1.second * r_2/r_1;

    diff.first = p2.first - x_t;
    diff.second = p2.second - y_t;
}

void MollerGEMSpatialRes::SetBeamEnergy(double e)
{
    beam_energy = e;
}

void MollerGEMSpatialRes::SetPoint1(pair<double, double> & p)
{
    p1 = p;
}

void MollerGEMSpatialRes::SetPoint2(pair<double, double> & p)
{
    p2 = p;
}

void MollerGEMSpatialRes::SetOrigin(pair<double, double> p)
{
    x_origin = p.first;
    y_origin = p.second;
}

void MollerGEMSpatialRes::SetGEM1Z(double z = 5320.)
{
    z_gem1 = z;
}

void MollerGEMSpatialRes::SetGEM2Z(double z = 5280.)
{
    z_gem2 = z;
}

pair<double, double> & MollerGEMSpatialRes::GetDiff()
{
    return diff;
}

double MollerGEMSpatialRes::EnergyFromAngle(double a)
{
    double e_elec = 0.511; // MeV
    //a = a/180. * PI;

    double _alpha = (beam_energy + e_elec) / TMath::Cos(a) / TMath::Sqrt(beam_energy * beam_energy - e_elec * e_elec);
    double alpha = _alpha * _alpha;

    return (alpha+1)/(alpha-1) * e_elec;
}

double MollerGEMSpatialRes::AngleFromEnergy(double e)
{
    // beam energy must in MeV
    double e_elec = 0.511; // MeV
    double numerator = e * beam_energy + e_elec * (e-beam_energy) - e_elec * e_elec;
    double denominator = TMath::Sqrt( (beam_energy * beam_energy - e_elec * e_elec) * (e * e - e_elec * e_elec));

    double cosangle = numerator / denominator;

    //return TMath::ACos(cosangle)*180./PI;
    return TMath::ACos(cosangle);
}

double MollerGEMSpatialRes::r( double x, double y)
{
    return TMath::Sqrt(x*x + y*y);
}
