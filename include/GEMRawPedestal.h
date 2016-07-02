#ifndef __GEMRAWPEDESTAL_H__
#define __GEMRAWPEDESTAL_H__

#include <TMath.h>
#include <TH1F.h>
#include "TCanvas.h"

#include "PRDMapping.h"

#include <iostream>
#include <cassert>

class GEMRawPedestal
{
public:
  GEMRawPedestal( map<int, map<int, vector<int> > > );
  ~GEMRawPedestal();

  void ClearMap();
  void ApvEventDecode();
  void ComputeApvPedestal(int, int);
  void ComputeEventPedestal();
  //int  IsChannelActive(int, int); //decoder already checked this, not necessary any more

  //void LoadPedestalData(const char* filename);
  //void SavePedestalRunHistos();
  void PrintEventPedestal(); // for debug

  Float_t GetStripNoise(int fecid, int adc_ch, int channelID);
  Float_t GetStripOffset(int fecid, int adc_ch, int channelID);
  
  Float_t GetStripMeanOffset(int fecid, int adc_ch, int channelID);

private:
  Float_t fApvHeaderLevel;
  Int_t NCH;
  Int_t fTimeSample;

  //Kondo's new fix
  Int_t fAPVID;
  TString fAPVStatus;
  //Int_t fCurrentEventTimeSample;
  Int_t fFECID;
  Int_t fADCCh;

  PRDMapping* mapping;

  map<int, map<int, vector<int> > > mSrsSingleEvent;
  //vector<int> vActiveAdcChannels;
  vector<int> vSingleApvData;
  map<int, map<int, vector<Float_t> > > mStripOffset;
  map<int, map<int, vector<Float_t> > > mStripRMS;

  vector<Float_t> vCommonModeOffset;
  vector<Float_t> vCommonModeOffset_split;
  multimap<Int_t, Float_t> mApvTimeBinData;


};

#endif
