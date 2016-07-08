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

#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TDirectory.h>

#include <stdio.h>

#include "GEMHistoManager.h"
#include "GEMConfigure.h"

class PRadReconstructor;
class PRadGEMReconstructor;
class PRadDSTParser;

class GEMPhysHandler : public GEMHistoManager
{
public:
  GEMPhysHandler();
  ~GEMPhysHandler();
  
  void ProcessAllFiles();
  int ProcessAllEvents(int evtID = -1);
  void SavePhysResults();

  void ProcessEp();
  void ProcessMoller();
  void CharactorizeGEM();
  void GeometryMollerRing(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct>&gem2);

private:
  string fileList[100];
  string filename;
  ifstream file;

  GEMConfigure config;

  PRadDataHandler *pHandler;
  PRadEvioParser *parser;
  PRadReconstructor *reconstruct;
  PRadGEMReconstructor *pGEMReconstructor;
  pRadDSTParser * pDSTParser;

  //temp test
  TH2F * hhTimeCorrelation;
  TH1F* hTimeDiff;

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
};

#endif
