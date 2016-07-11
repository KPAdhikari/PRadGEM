#ifndef __GEMDATASTRUCT__H__
#define __GEMDATASTRUCT__H__

struct GEMClusterStruct
{
    float x;
    float y;
    float x_charge;
    float y_charge;
    float energy;
    float z;

    GEMClusterStruct(float xi, float yi, 
               float cix = 0., float ciy = 0., 
	       float ei = 0) 
    : x(xi), y(yi), x_charge(cix), y_charge(ciy), energy(ei) {}

    void SetEnergy(float e) {energy = e;}
    void SetX(float xp) {x = xp;}
    void SetY(float xp) {y = xp;}
    void SetXCharge(float xp) {x_charge = xp;}
    void SetYCharge(float xp) {y_charge = xp;}
    void SetZ(float xp) { z = xp;}

};

#endif
