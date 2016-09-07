#include "PRadEP.h"
#include <TMath.h>
#include <iostream>

#define UNDEFINED_VALUE -1000.
#define PI 3.1415926

using namespace std;

PRadEP::PRadEP()
{
    beam_energy = 2100;
}

PRadEP::~PRadEP()
{
}

void PRadEP::SetData(vector<GEMClusterStruct> *fgem)
{
    gem = fgem;
}

void PRadEP::SetBeamEnergy(double & e)
{
    beam_energy = e;
}

void PRadEP::SetBeamEnergy(double && e)
{
    beam_energy = e;
}

bool PRadEP::PassCut()
{
    // preliminary cut
    if(gem->size() != 1)
        return false;
    else if( ! EnergyCut( gem->at(0).energy ) )
        return false;
    else 
        return true;
}

void PRadEP::ProcessData(vector<GEMClusterStruct> *fgem)
{
    Reset();
    SetData(fgem);
    if( ! PassCut() )
        return;
    Process();
}

void PRadEP::Process()
{
    double theta = TMath::ATan( r(gem->at(0).x, gem->at(0).y)/gem->at(0).z); 
    scatt_angle = RadToDec(theta);
    energy_angle.emplace_back(gem->at(0).energy, scatt_angle);
    positions.emplace_back(gem->at(0).x, gem->at(0).y);
    chamber_id = gem->at(0).chamber_id;

    // q square
    double numerator = 4.0*(beam_energy/1000)*(beam_energy/1000)*TMath::Sin(theta/2)*TMath::Sin(theta/2);
    double denominator = 1+ (2*(beam_energy/1000)/0.938)*TMath::Sin(theta/2)*TMath::Sin(theta/2);
    q_square = numerator/denominator;
}

void PRadEP::Reset()
{
    q_square = UNDEFINED_VALUE;
    scatt_angle = UNDEFINED_VALUE;
    energy_angle.clear();
    positions.clear();
    hycal_match.clear();
}

bool PRadEP::EnergyCut( double && e)
{
    //return ( (e>beam_energy*0.85) && (e<beam_energy*1.15) );
    return true;
}

double PRadEP::r(double & x, double & y)
{
    return TMath::Sqrt(x*x + y*y);
}

double PRadEP::r(double && x, double && y)
{
    return TMath::Sqrt(x*x + y*y);
}

double PRadEP::RadToDec(double && v)
{
    return v*180/PI;
}

double PRadEP::RadToDec(double & v)
{
    return v*180/PI;
}

double & PRadEP::ScattAngle()
{
    return scatt_angle;
}

double & PRadEP::QSquare()
{
    return q_square;
}

vector<pair<double, double> > & PRadEP::EnergyAngle()
{
    return energy_angle;
}

vector<pair<double, double> > & PRadEP::Positions()
{
    return positions;
}

void PRadEP::SetEvtID(unsigned int id)
{
    evt_id = id;
}

unsigned int PRadEP::GetEvtID()
{
    return evt_id;
}

int & PRadEP::GetChamberID()
{
    return chamber_id;
}

void PRadEP::SetChamberID(int id)
{
    chamber_id = id;
}

void PRadEP::SetHyCalMatch(vector<HyCalHit> hycal)
{
    hycal_match = hycal;
}

vector<HyCalHit> & PRadEP::GetHyCalMatch()
{
    return hycal_match;
}
