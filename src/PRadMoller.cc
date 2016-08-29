//==============================================================================//
// Note: This class only handles one moller event per trigger,                  //
//       adapted for PRad, which has a low event rate.                          //
//                                                                              //
// Xinzhan Bai                                                                  //
// 08/09/2016                                                                   //
//==============================================================================//
#include "PRadMoller.h"
#include "MollerGEMSpatialRes.h"
#include <TMath.h>
#include <iostream>
#include <cassert>

#define PI 3.1415926
#define Undefined_Value -1000.

using namespace std;

PRadMoller::PRadMoller()
{
    beam_energy = 1100;//1.1GeV
    previous_positions.reserve(2);
    previous_positions.emplace_back(Undefined_Value, Undefined_Value);
    previous_positions.emplace_back(Undefined_Value, Undefined_Value);
    gem_pos_res = new MollerGEMSpatialRes();
}

PRadMoller::~PRadMoller()
{
}

void PRadMoller::SetBeamEnergy(double &e)
{
    beam_energy = e;
    gem_pos_res->SetBeamEnergy(e);
}

void PRadMoller::SetBeamEnergy(double &&e)
{
    beam_energy = e;
    gem_pos_res -> SetBeamEnergy(e);
}

void PRadMoller::SetData(vector<GEMClusterStruct> *fgem)
{
    gem = fgem;
}

void PRadMoller::Reset()
{
    scatt_angle1 = Undefined_Value;
    scatt_angle2 = Undefined_Value;
    open_angle = Undefined_Value;
    coplanarity = Undefined_Value;
    moller_center.first = Undefined_Value;
    moller_center.second = Undefined_Value;
    energy_angle.clear();
    positions.clear();
}

bool PRadMoller::PassCut()
{
    // preliminary cut
    if(gem->size() != 2)
        return false;
    if( ! EnergyCut(gem->at(0).energy,  gem->at(1).energy ) )
	return false;
    //if( ! QuandrantsCut() )
    //	return false;
    return true;
}

bool PRadMoller::QuandrantsCut()
{
    //--------------------------------
    // intent for this cut:
    //   if the two clusters are in 
    //   opposite quadrants, then this
    //   is a moller event, otherwise 
    //   it is not a moller event
    //
    //  is this cut really necessary?
    //---------------------------------
    if( (gem->at(0).x * gem->at(1).x <= 0) &&
	(gem->at(0).y * gem->at(1).y <= 0) )
	return true;
    else 
        return false;
}

bool PRadMoller::EnergyCut(double && e1, double &&e2)
{
    // 3 sigma(3%) cut
    return ((e1+e2) > beam_energy *0.91) && ( (e1+e2)< beam_energy*1.09);
    //return true;
}

void PRadMoller::ProcessData(vector<GEMClusterStruct> *fgem)
{
    Reset();
    SetData(fgem);
    if( ! PassCut() )
	return;
    Process();
}

void PRadMoller::Process()
{
    scatt_angle1 = RadToDec( TMath::ATan(r(gem->at(0).x, gem->at(0).y)/gem->at(0).z) );
    scatt_angle2 = RadToDec( TMath::ATan(r(gem->at(1).x, gem->at(1).y)/gem->at(1).z) );
    open_angle = scatt_angle1 + scatt_angle2;
    coplanarity = GetCoplanarity();    
    energy_angle.emplace_back(gem->at(0).energy, scatt_angle1);
    energy_angle.emplace_back(gem->at(1).energy, scatt_angle2);
    positions.emplace_back(gem->at(0).x, gem->at(0).y);
    positions.emplace_back(gem->at(1).x, gem->at(1).y);
    GetMollerCenter();
    //-----------------------------------------------------------
    // computed moller center history:
    //
    // date before 08/16/2016: 1.473, 0.601
    // 08/17/2016: 1.448, 1.289
    //-----------------------------------------------------------
    gem_pos_res->SetOrigin(make_pair<double, double>(1.448, 1.289));
    gem_pos_res->Process(positions);
}

double PRadMoller::GetCoplanarity()
{
    double a1 = GetXYSlopeAngle(gem->at(0).x, gem->at(0).y);
    double a2 = GetXYSlopeAngle(gem->at(1).x, gem->at(1).y);
    if( a1 > a2 )
	return a1-a2-180.;
    else if(a1 < a2)
	return a2-a1-180.;
    else{
	cout<<"Error: two hits on hycal matching one same hit on gem."
	    <<endl;
    }
}

MollerGEMSpatialRes * PRadMoller::GetSpatialResHandler()
{
    return gem_pos_res;
}

double PRadMoller::GetXYSlopeAngle(double && x, double && y)
{
    int quadrant = GetQuandrant(x, y);
    double slope = 0.;
    switch( quadrant)
    {
	case 1: {
		    slope = TMath::Abs(y/x);
		    slope = TMath::ATan(slope);
		    break;
		}
	case 2: {
		    slope = TMath::Abs(y/x);
		    slope = TMath::ATan(slope);
		    slope = PI - slope;
		    break;
		}
	case 3: {
		    slope = TMath::Abs(y/x);
		    slope = TMath::ATan(slope);
		    slope = PI + slope;
		    break;
		}
	case 4: {
		    slope = TMath::Abs(y/x);
		    slope = TMath::ATan(slope);
		    slope = 2*PI - slope;
		    break;
		}
	case 0: {
		    if( x == 0) {
			if( y>0)
			    slope = PI/2;
			else
			    slope = 3*PI/2;
		    }
		    else if(y == 0) {
			if(x>0)
			    slope = 0;
			else
			    slope = PI;
		    }
		    else {
			cout<<"Get Slope Error..."<<endl;
		    }
		    break;
		}
	default: {
		     slope = Undefined_Value;
		     cout<<" GetQuandrantError..."<<endl;
		     break;
		 }
    };
    return slope*180/PI;
}

double PRadMoller::r(double & x, double & y)
{
    return TMath::Sqrt(x*x + y*y);
}

double PRadMoller::r(double && x, double && y)
{
    return TMath::Sqrt(x*x + y*y);
}

double PRadMoller::RadToDec(double && v)
{
    return v*180/PI;
}

double PRadMoller::RadToDec(double & v)
{
    return v*180/PI;
}

int PRadMoller::GetQuandrant(double & x, double & y)
{
    if( (x>0) && (y>0))
	return 1;
    else if( (x<0) && (y>0))
	return 2;
    else if( (x<0) && (y<0) )
	return 3;
    else if( (x>0) && (y<0) )
	return 4;
    else
	return 0;
    /* point lying on axis */
}

int PRadMoller::GetQuandrant(double && x, double && y)
{
    if( (x>0) && (y>0))
	return 1;
    else if( (x<0) && (y>0))
	return 2;
    else if( (x<0) && (y<0) )
	return 3;
    else if( (x>0) && (y<0) )
	return 4;
    else
	return 0;
    /* point lying on axis */
}

double & PRadMoller::ScattAngle1()
{
    return scatt_angle1;
}

double & PRadMoller::ScattAngle2()
{
    return scatt_angle2;
}

double & PRadMoller::OpenAngle()
{
    return open_angle;
}

double & PRadMoller::Coplanarity()
{
    return coplanarity;
}

vector<pair<double, double> > & PRadMoller::EnergyAngle()
{
    return energy_angle;
}

vector<pair<double, double> > & PRadMoller::Positions()
{
    return positions;
}

pair<double, double> & PRadMoller::MollerCenter()
{
    return moller_center;
}

void PRadMoller::GetMollerCenter()
{
    assert(previous_positions.size() == 2);
    if(positions.size() != 2)
	return;
    if(previous_positions[0].first == Undefined_Value){
	previous_positions.clear();
	previous_positions.push_back(positions[0]);
	previous_positions.push_back(positions[1]);
	return;
    }
    else {
	GetIntersection();
	previous_positions.clear();
	previous_positions.push_back(positions[0]);
	previous_positions.push_back(positions[1]);
    }
}

void PRadMoller::GetIntersection()
{
    double _x1 = previous_positions[0].first;
    double _y1 = previous_positions[0].second;
    double _x2 = previous_positions[1].first;
    double _y2 = previous_positions[1].second;
    double x1 = positions[0].first;
    double y1 = positions[0].second;
    double x2 = positions[1].first;
    double y2 = positions[1].second;

    double a1 = _y1 - _y2;
    double b1 = _x2 - _x1;
    double c1 = _y2*_x1 - _y1*_x2;

    double a2 = y1 - y2;
    double b2 = x2 - x1;
    double c2 = y2*x1 - y1*x2;

    moller_center.second = (a1*c2 - a2*c1 )/(a2*b1 - a1*b2);
    moller_center.first = (b1*c2 - b2*c1)/(a1*b2 - a2*b1);

    //debug
    //cout<<"("<<_x1<<","<<_y1<<") ("<<_x2<<", "<<_y2<<")"<<endl;
    //cout<<"("<<x1<<","<<y1<<") ("<<x2<<", "<<y2<<")"<<endl;
    //cout<<"("<<moller_center.first<<", "<<moller_center.second<<")"<<endl;
}
