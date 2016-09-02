#ifndef PRAD_EP_H
#define PRAD_EP_H
#include <vector>
#include "GEMDataStruct.h"

class PRadEP
{
public:
    PRadEP();
    ~PRadEP();
    void SetData(std::vector<GEMClusterStruct> * fgem);
    void SetBeamEnergy(double & );
    void SetBeamEnergy(double &&);
    bool PassCut();
    bool EnergyCut(double &&);
    void ProcessData(std::vector<GEMClusterStruct> *fgem);
    void Reset();
    void Process();
    double r(double &, double &);
    double r(double &&, double &&);
    double RadToDec(double &);
    double RadToDec(double &&);
    // results
    double & ScattAngle();
    double & QSquare();
    std::vector<std::pair<double, double> > & EnergyAngle();
    std::vector<std::pair<double, double> > & Positions();

    void SetEvtID( unsigned int );
    unsigned int GetEvtID();
    int & GetChamberID();
    void SetChamberID(int );

private:
    unsigned int evt_id;
    std::vector<GEMClusterStruct> *gem;
    double beam_energy;

    double scatt_angle;
    double q_square;
    std::vector< std::pair<double, double > > energy_angle;
    std::vector< std::pair<double, double > > positions;
    int chamber_id;
};

#endif
