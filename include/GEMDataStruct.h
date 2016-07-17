#ifndef __GEMDATASTRUCT__H__
#define __GEMDATASTRUCT__H__

#include <TObject.h>

class GEMClusterStruct : public TObject
{
public:
    float x;
    float y;
    float x_charge;
    float y_charge;
    float energy;
    float z;
    int x_size;
    int y_size;

    GEMClusterStruct();
    ~GEMClusterStruct();

    GEMClusterStruct(float xi, float yi, 
               float cix = 0., float ciy = 0., 
	       float ei = 0) 
    : x(xi), y(yi), x_charge(cix), y_charge(ciy), energy(ei) 
    {
      z = 0;
      x_size = 0;
      y_size = 0;
    }

    GEMClusterStruct( const GEMClusterStruct & gem_cluster)
    : x(gem_cluster.x), y(gem_cluster.y),  
      x_charge(gem_cluster.x_charge), 
      y_charge(gem_cluster.y_charge),
      energy(gem_cluster.energy), 
      z(gem_cluster.z), 
      x_size(gem_cluster.x_size),
      y_size(gem_cluster.y_size){}

    void SetEnergy(float e) {energy = e;}
    void SetX(float xp) {x = xp;}
    void SetY(float xp) {y = xp;}
    void SetXCharge(float xp) {x_charge = xp;}
    void SetYCharge(float xp) {y_charge = xp;}
    void SetZ(float xp) { z = xp;}

ClassDef(GEMClusterStruct, 1)
};

class HyCalClusterStruct : public TObject
{
public:
    float x;
    float y;
    float energy;

    HyCalClusterStruct();
    ~HyCalClusterStruct();

    HyCalClusterStruct(float xi, float yi, float energyi)
        : x(xi), y(yi), energy(energyi) {}

    HyCalClusterStruct( const HyCalClusterStruct & hycal)
        : x(hycal.x), y(hycal.y), energy(hycal.energy) {}

ClassDef(HyCalClusterStruct, 1)

};

#endif
