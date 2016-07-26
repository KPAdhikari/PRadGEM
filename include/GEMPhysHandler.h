#ifndef __GEMPHYS_HANDLER_H__
#define __GEMPHYS_HANDLER_H__

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <string>
#include <cassert>
#include <algorithm>
#include <unordered_map>

#include <arpa/inet.h>

#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TDirectory.h>

#include <stdio.h>

//evio
#include "evioUtil.hxx"
#include "evioFileChannel.hxx"

#include "GEMRawDecoder.h"
#include "PRDMapping.h"
#include "GEMRawPedestal.h"
#include "GEMHistoManager.h"
#include "GEMConfigure.h"
#include "GEMPedestal.h"
#include "GEMDataStruct.h"

//decoding data after zero suppression
#include "GEMZeroHitDecoder.h"

// chao
#include "PRadDataHandler.h"
#include "PRadEvioParser.h"

class PRadGEMTree;
class PRadReconstructor;

class GEMPhysHandler : public GEMHistoManager
{
public:
  GEMPhysHandler();
  ~GEMPhysHandler();
  
  void ProcessAllFiles();
  int ProcessAllEvents(int evtID = -1);
  void SavePhysResults();

  template<class T> void ProcessEp(T *hit_decoder);
  template<class T> void ProcessMoller(T *hit_decoder);
  template<class T> void CharactorizeGEM(T *hit_decoder);
  template<class T> void CharactorizeHyCal(T *hit_decoder);
  template<class T> void ComputeGEMOffsets(T *hit_decoder);
  template<class T> void GetGEMClusterMapHyCalCoor(T *hit_decoder);

  void GeometryMollerRing(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct>&gem2);

  template<class T> void ProcessMollerAfterCorrection(T *hit_decoder);
  template<class T> void EvalMatchMech(T *hit_decoder);

  int GEMHyCalPosMatch(int i, vector<GEMClusterStruct> &gem, vector<HyCalHit> *pHHit);
  int HyCalGEMPosMatch( vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2, vector<HyCalHit> *pHHit);
  void GetCutFlags(int nth_event, vector<uint32_t> * vec, int & scin_flag, int & hycal_flag);

  //init
  void InitSetup();
  void InitPedestal();
  void InitVariables();
  void InitHandler();
  void InitIntermediateMollerCenterVariable();
  void BookTimingHistos();
  
  void InitTDCGroup();
  int GetTDCGroup(const string &);

  // moller energy vs scattering angle formula
  double MollerAngleFromEnergy(double E);
  double MollerEnergyFromAngle(double theta);

  // sectorize gem
  int GetSectorIndex( double x, double y);
  void WriteSectorEff();

private:
  vector<int> vSRSSingleEventData;
  vector<int> vSRSZeroEventData;
  unordered_map<string, int> TDC_Map; 
  set<int> FECs;  // FEC ID

  string fileList[100];
  string filename;
  ifstream file;

  int nElectron;  
  // tdc 126 or 127 cut, converted events number
  int nElectron_126;
  int nElectron_127;
  int neff;
  int neff_after_match;
 
  int nScinEvents;
  int nHyCalEvents;
  // totoal number of events 
  unsigned long nTotalEvents;

  // calculte gem efficiency using ep 
  // and moller events from production
  // runs, during the matching process
  double GEMMollerElectronQuantity;
  double HyCalMollerElectronQuantity;
  double GEMEpElectronQuantity;
  double HyCalEpElectronQuantity;
  // sectorize GEM
  /*
  float x_sector[14] ={-826.5, -345.9, -341.8, -162.5, -159.8, -37.7,
                        39.1,  159.9,  163.9,  343.2,  345.9, 403.0,
		       404.4, 826.5}; 
  float y_sector[16] ={-885.0, -410.9, -408.1, -205.5, -202.7, -39.0,
                        -1.3,    1.3,   39.0,  202.7,  205.9, 302.7,
		       312.0,  408.2,  410.8, 885.0};
  */
  float x_sector[14] ={-826.5, -355.9, -331.8, -172.5, -149.8, -47.7,
                        49.1,  149.9,  173.9,  333.2,  355.9, 393.0,
		       414.4, 826.5}; 
  float y_sector[16] ={-885.0, -420.9, -398.1, -215.5, -192.7, -49.0,
                        -11.3,    11.3,   49.0,  192.7,  215.9, 292.7,
		       322.0,  398.2,  420.8, 885.0};

  double gem_moller_quantity[70];
  double hycal_moller_quantity[70];
  double gem_ep_quantity[70];
  double hycal_ep_quantity[70];
  // circular sector
  float r_sector[50];
  double gem_r_sec_ep_quantity[50];
  double hycal_r_sec_ep_quantity[50];

  GEMRawDecoder *fRawDecoder;
  GEMConfigure *config;
  GEMPedestal *ped;
  // chao
  PRadDataHandler *pHandler;
  PRadEvioParser *parser;
  PRadReconstructor *reconstruct;

  double timing_test;
  vector<HyCalHit> * pHyCalHit;  
  double totalEnergyDeposit;
  double beamEnergy;
  double beamEnergyCut;

  //compute intersection points between two moller events
  float px1, py1, px2, py2;
  float cx1, cy1, cx2, cy2;
  float px1_c, py1_c, px2_c, py2_c;
  float cx1_c, cy1_c, cx2_c, cy2_c;

  // coordinate transfer
  // origin moved from chamber center to beam hole area
  /* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   * overlapping area: 44mm
   * x side length: 550.4
   * overlapping area starting with: 550.4/2 -44 = 231.2
   *
   * origin trasfer to beam hole center:
   *    transfered distance: 550.4/2 - 44/2 = 253.2
   *    GEM1 coordinate transfer: x1 = x1 - 253.2; y1 = y1
   *    GEM2 coordinate transfer: x2 = 253.2 - x2; y2 = -y2
   *xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  float O_Transfer;
  float OverlapStart;

  float Z_gem1;
  float Z_gem2;
  float Z_hycal;
  float Delta;

  // save results
  fstream outfile;
  PRadGEMTree *rst_tree;

};

#endif
