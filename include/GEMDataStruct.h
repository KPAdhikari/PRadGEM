#ifndef __GEMDATASTRUCT__H__
#define __GEMDATASTRUCT__H__

struct EventHeader
{
    unsigned int length;
    unsigned char num;
    unsigned char type;
    unsigned short tag;
};

struct GEMClusterStruct
{
    float x;
    float y;
    float x_charge;
    float y_charge;
    float energy;
    float z;
    int x_size;
    int y_size;
    int chamber_id;

    GEMClusterStruct(){};
    ~GEMClusterStruct(){};

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
    void SetXClusterSize( int s) {x_size = s;}
    void SetYClusterSize( int s) {y_size = s;}
};

struct HyCalClusterStruct
{
    float x;
    float y;
    float energy;

    HyCalClusterStruct(){};
    ~HyCalClusterStruct(){};

    HyCalClusterStruct(float xi, float yi, float energyi)
	: x(xi), y(yi), energy(energyi) {}

    HyCalClusterStruct( const HyCalClusterStruct & hycal)
	: x(hycal.x), y(hycal.y), energy(hycal.energy) {}
};

struct EpicsStruct
{
    double MFC_Flow;
    double Cell_Gas_T;
    double Tank_P_P;
    double Chamber_P;
    double Cell_P;
    double Cell_Body_T;
    double hallb_ptrans_y2_encoder;
    double hallb_ptrans_y1_encoder;
    double hallb_ptrans_x_encoder;
    double ptrans_y;
    double ptrans_x;
    double hallb_IPM2H01_CUR;
    double hallb_IPM2H01_YPOS;
    double hallb_IPM2H01_XPOS;
    double hallb_IPM2C24A_CUR;
    double hallb_IPM2C24A_XPOS;
    double hallb_IPM2C21A_CUR;
    double hallb_IPM2C21A_YPOS;
    double hallb_IPM2C21A_XPOS;
    double hallb_IPM2C24A_YPOS;
    double scaler_calc1;
    double VCG2H02A;
    double VCG2H01A;
    double VCG2H00A;
    double VCG2C24A;
    double VCG2C21A;
    double VCG2C21;
    double MBSY2C_energy;

    EpicsStruct()
	: MFC_Flow(-1000), Cell_Gas_T(-1000), Tank_P_P(-1000),
	Chamber_P(-1000), Cell_P(-1000), Cell_Body_T(-1000),
	hallb_ptrans_y2_encoder(-1000), hallb_ptrans_y1_encoder(-1000),
	hallb_ptrans_x_encoder(-1000), ptrans_y(-1000), ptrans_x(-1000),
	hallb_IPM2H01_CUR(-1000), hallb_IPM2H01_YPOS(-1000),
	hallb_IPM2H01_XPOS(-1000), hallb_IPM2C24A_CUR(-1000),
	hallb_IPM2C24A_XPOS(-1000), hallb_IPM2C21A_CUR(-1000),
	hallb_IPM2C21A_YPOS(-1000), hallb_IPM2C21A_XPOS(-1000),
	hallb_IPM2C24A_YPOS(-1000), scaler_calc1(-1000),
	VCG2H02A(-1000), VCG2H01A(-1000), VCG2H00A(-1000),
	VCG2C24A(-1000), VCG2C21A(-1000), VCG2C21(-1000),
	MBSY2C_energy(-1000)
    {};

    void Clear()
    {
	MFC_Flow = -1000.; Cell_Gas_T = -1000.; Tank_P_P = -1000.;
	Chamber_P = -1000.; Cell_P = -1000.;  Cell_Body_T = -1000.;
	hallb_ptrans_y2_encoder = -1000.; hallb_ptrans_y1_encoder = -1000.;
	hallb_ptrans_x_encoder = -1000.; ptrans_y = -1000.; ptrans_x = -1000.;
	hallb_IPM2H01_CUR = -1000.; hallb_IPM2H01_YPOS = -1000.;
	hallb_IPM2H01_XPOS = -1000.; hallb_IPM2C24A_CUR = -1000.;
	hallb_IPM2C24A_XPOS = -1000.; hallb_IPM2C21A_CUR = -1000.;
	hallb_IPM2C21A_YPOS = -1000.; hallb_IPM2C21A_XPOS = -1000.;
	hallb_IPM2C24A_YPOS=-1000.; scaler_calc1=-1000.;
	VCG2H02A = -1000.; VCG2H01A = -1000.; VCG2H00A = -1000.;
	VCG2C24A = -1000.; VCG2C21A =-1000.; VCG2C21 =-1000.;
	MBSY2C_energy = -1000.;
    }
};

#endif
