#ifndef __GEMPEDESTAL_H__
#define __GEMPEDESTAL_H__

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <string>
#include <cassert>

#include <arpa/inet.h>

#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>

#include <stdio.h> //for getchar()

//evio
#include "evioUtil.hxx"
#include "evioFileChannel.hxx"

#include "GEMRawDecoder.h"
#include "PRDMapping.h"
#include "GEMRawPedestal.h"

#include "GEMConfigure.h"

//using namespace std;
//using namespace evio;

class GEMPedestal
{
public:
  GEMPedestal(string str);
  ~GEMPedestal();

  void BookHistos();

  int ProcessAllEvents(int evtID = -1);
  void ComputePedestal();
  void SavePedestalFile();
  void LoadPedestal();
  void Delete();
  vector<Float_t> GetAPVNoises(Int_t);
  vector<Float_t> GetAPVOffsets(Int_t);

private:
  vector<int> vSRSSingleEventData;
  map<int, map<int, vector<int> > > mAPVRawTSs;

  set<int> FECs;  // FEC ID
  Int_t NCH;
  Int_t nNbofAPVs;

  string filename;
  string pedestal_file;

  ifstream file;

  vector<TH1F*> vStripOffsetHistos;
  vector<TH1F*> vStripNoiseHistos;
  vector<TH1F*> vApvPedestalOffset;
  vector<TH1F*> vApvPedestalNoise;

  TH1F* hAllStripNoise;
  TH1F* hAllXStripNoise;
  TH1F* hAllYStripNoise;

  // not care for now, will implement later if interested
  //TH1F* hAllStripOffset;
  //TH1F* hAllXStripOffset;
  //TH1F* hAllYStripOffset;

  PRDMapping * mapping;

  //GEMRawDecoder *fRawDecoder;
  GEMConfigure config;
};

#endif
