#ifndef __GEMINPUT_HANDLER_H__
#define __GEMINPUT_HANDLER_H__

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
#include "GEMPedestal.h"
#include "GEMConfigure.h"
#include "GEMDataStruct.h"

//using namespace std;
//using namespace evio;

class GEMInputHandler
{
public:
  GEMInputHandler(string);
  ~GEMInputHandler();

  int ProcessAllEvents(int evtID = -1);
  //bool FECEventChooser(const evioDOMNodeP pNode);

private:
  vector<int> vSRSSingleEventData;
  vector<int> vROCZeroSupData;
  map<int, map<int, TH1F*> > mAPVRawHistos;
  map<int, map<int, vector<int> > > mAPVRawTSs;

  set<int> FECs;  // FEC ID

  string filename;
  ifstream file;

  GEMPedestal *ped;
  GEMConfigure configure;

  GEMRawDecoder *fRawDecoder;
};

#endif
