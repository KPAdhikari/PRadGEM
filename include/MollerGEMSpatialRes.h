#ifndef MOLLER_GEM_SPATIAL_RES_H
#define MOLLER_GEM_SPATIAL_RES_H
#include <vector>
#include "GEMDataStruct.h"

class MollerGEMSpatialRes
{
public:
    MollerGEMSpatialRes();
    ~MollerGEMSpatialRes();

    double EnergyFromAngle(double);
    double AngleFromEnergy(double);
    void SetBeamEnergy(double);
    void SetPoint1(std::pair<double, double>&);
    void SetPoint2(std::pair<double, double>&);
    void SetOrigin(std::pair<double, double>);
    void SetZ(double);
    double r(double, double);
    void Reset();
    void Process(std::vector<std::pair<double, double> > );
    void _Process();
    std::pair<double, double>& GetDiff();

private:
    double beam_energy;
    double x_origin;
    double y_origin;
    double z_gem;
    std::pair<double, double> p1;
    std::pair<double, double> p2;
    std::pair<double, double> diff;
};

#endif
