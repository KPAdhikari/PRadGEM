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
  
  void InitTDCGroup();
  int GetTDCGroup(const string &);

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
  int nTotalEvents;

  // calculte gem efficiency using ep 
  // and moller events from production
  // runs, during the matching process
  double GEMMollerElectronQuantity;
  double HyCalMollerElectronQuantity;
  double GEMEpElectronQuantity;
  double HyCalEpElectronQuantity;

  GEMRawDecoder *fRawDecoder;
  GEMConfigure config;
  GEMPedestal *ped;
  // chao
  PRadDataHandler *pHandler;
  PRadEvioParser *parser;
  PRadReconstructor *reconstruct;

  //temp test
  TH2F * hhTimeCorrelation;
  TH1F* hTimeDiff;
  double timing_test;
  vector<HyCalHit> * pHyCalHit;  
  //vector<HyCalHit>  HyCalHit; 
  double totalEnergyDeposit;
  double beamEnergy;

  //compute intersection points between two moller events
  float px1, py1, px2, py2;
  float cx1, cy1, cx2, cy2;
  float px1_c, py1_c, px2_c, py2_c;
  float cx1_c, cy1_c, cx2_c, cy2_c;

  // coordinate transfer
  // origin moved from chamber center to beam hole area
  /*
   * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   * overlapping area: 44mm
   * x side length: 550.4
   * overlapping area starting with: 550.4/2 -44 = 231.2
   *
   * origin trasfer to beam hole center:
   *    transfered distance: 550.4/2 - 44/2 = 253.2
   *    GEM1 coordinate transfer: x1 = x1 - 253.2; y1 = y1
   *    GEM2 coordinate transfer: x2 = 253.2 - x2; y2 = -y2
   * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   */

  float O_Transfer;
  float OverlapStart;

  // save results
  fstream outfile;

};

#endif
