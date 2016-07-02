#ifndef _GEMONLINEHITDECODER_H__
#define _GEMONLINEHITDECODER_H__

#include "GEMRawDecoder.h"
#include "PRDMapping.h"
#include "GEMHit.h"
#include "GEMPedestal.h"
#include "GEMCluster.h"
#include "GEMDataStruct.h"

#include <stdint.h>

#include "TString.h"
#include "TObject.h"
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include "TH1F.h"
#include "TList.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMath.h"

#include <map>
#include <list>
#include <vector>
#include <numeric>

class GEMOnlineHitDecoder
{
public:
  GEMOnlineHitDecoder(uint32_t *buf, int Size, GEMPedestal *pPed);
  ~GEMOnlineHitDecoder();

  void ProcessEvent();
  void EventHandler(map<int, map<int, vector<int> > > , GEMPedestal*);

  map< TString, list<GEMHit*> > GetListOfHitsFromPlanes();
  TH1F* GetHit(TString plane);

  map< TString, list<GEMHit*> > GetListOfHitsCleanFromPlanes();
  TH1F* GetCleanHit(TString plane);

  Bool_t IsADCchannelActive() ;

  void APVEventDecoder();
  void APVEventSplitChannelsDecoder();

  // compute cluster information 
  void ComputeClusters();
  void DeleteClustersInPlaneMap() ;
  map < TString, list <GEMCluster * > > GetListOfClustersFromPlanes() { return  fListOfClustersCleanFromPlane;  }

  TH1F* GetCluster(TString str);

  void GetClusterGEM(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2);
  void GetClusterHyCalCutMode(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2);
  void GetClusterHyCalPlusMode(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2);
  void GetClusterBeamLine(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2);

  void FillHistos(TH1F**, TH1F**, TH1F**, TH1F**);

private:
  uint32_t *buf;
  int fSize;
  Int_t NCH;

  int fFECID;
  int fADCChannel;
  int fAPVID;
  int fAPVHeaderLevel;
  int fAPVKey;
  TString fIsHitMaxOrTotalADCs;
  vector<int> fActiveADCchannels;
  PRDMapping * fMapping;
  GEMPedestal * ped;
  set<int> FECs;
  TString fAPVStatus;

  map<Int_t, GEMHit*> fListOfHits;
  map<Int_t, GEMHit*> fListOfHitsClean;
  map<TString, list<GEMHit*> > fListOfHitsFromPlane;
  map<TString, list<GEMHit*> > fListOfHitsCleanFromPlane;

  string pedestal_file;
  int fZeroSupCut;
  int fSaveRawHit;
  int fStartData;
  int fTimeSample;

  map<int, map<int, vector<int> > > mSrsSingleEvent;
  vector<int> fRawData16bits;
  vector <Float_t> fPedestalNoises, fPedestalOffsets;
  vector <Float_t> fPedestalNoises_1stSet, fPedestalNoises_2ndSet, fPedestalOffsets_1stSet, fPedestalOffsets_2ndSet;

  TString fIsClusterMaxOrTotalADCs;
  Bool_t fIsGoodClusterEvent;
  Int_t fMinClusterSize, fMaxClusterSize;

  map<TString, list<GEMCluster*> > fListOfClustersCleanFromPlane;
};


#endif
