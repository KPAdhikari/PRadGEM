#ifndef PRADMOLLER_H
#define PRADMOLLER_H

#include "GEMDataStruct.h"
#include <vector>

class MollerGEMSpatialRes;

/* Note: This class only handles one moller event per trigger */

class PRadMoller
{
public:
    PRadMoller();
    ~PRadMoller();
    void SetData(std::vector<GEMClusterStruct> * fgem);
    void SetBeamEnergy(double & );
    void SetBeamEnergy(double &&);
    bool PassCut();
    bool EnergyCut(double &&, double &&);
    bool QuandrantsCut();
    void ProcessData(std::vector<GEMClusterStruct> *fgem);
    void Reset();
    void Process();
    double r(double &, double &);
    double r(double &&, double &&);
    double RadToDec(double &);
    double RadToDec(double &&);
    int GetQuandrant(double &&, double &&);
    int GetQuandrant(double &, double &);
    double GetCoplanarity();
    double GetXYSlopeAngle(double &&, double &&);
    void GetMollerCenter();
    void GetIntersection();
    // following get results
    double & ScattAngle1();
    double & ScattAngle2();
    double & OpenAngle();
    double & Coplanarity();
    std::vector<std::pair<double, double> > & EnergyAngle();
    std::vector<std::pair<double, double> > & Positions();
    std::pair<double, double> & MollerCenter();
    MollerGEMSpatialRes * GetSpatialResHandler();

private:
    MollerGEMSpatialRes * gem_pos_res;
    std::vector<GEMClusterStruct> *gem;
    double beam_energy;

    double scatt_angle1;
    double scatt_angle2;
    double open_angle;
    double coplanarity;
    std::vector< std::pair<double, double > > energy_angle;
    std::vector< std::pair<double, double > > positions;
    // moller center
    // note: moller center is on z=5230 plane
    std::pair<double , double > moller_center;
    std::vector< std::pair<double, double> > previous_positions;
};

#endif
