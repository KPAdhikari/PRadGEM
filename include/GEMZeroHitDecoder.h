#ifndef _GEMZEROHITDECODER_H__
#define _GEMZEROHITDECODER_H__

#include "PRDMapping.h"
#include "GEMHit.h"
#include "GEMCluster.h"
#include "GEMPedestal.h"
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


class GEMZeroHitDecoder
{
public:
  GEMZeroHitDecoder(uint32_t *buf, int Size, GEMPedestal *pPed);
  GEMZeroHitDecoder(uint32_t *buf, int Size);
  ~GEMZeroHitDecoder();

  void ProcessEvent();
  void EventHandler();

  map< TString, list<GEMHit*> > GetListOfHitsZeroFromPlanes();
  TH1F* GetZeroHit(TString plane);

  void Cut();

  // compute cluster information 
  void ComputeClusters();
  void ComputeClustersFineTune();
  void DeleteClustersInPlaneMap() ;
  map < TString, list <GEMCluster * > > GetListOfClustersFromPlanes() { 
      return  fListOfClustersZeroFromPlane;  
  }
 
  TH1F* GetCluster(TString str);

  void GetClusterGEM(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2);
  void GetClusterHyCalCutMode(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2);
  void GetClusterHyCalPlusMode(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2);
  void GetClusterBeamLine(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2);

  void FillHistos(TH1F**, TH1F**, TH1F**, TH1F**);

private:
  TString fIsClusterMaxOrTotalADCs;
  Bool_t fIsGoodClusterEvent;
  Int_t fMinClusterSize, fMaxClusterSize;

  map<TString, list<GEMCluster*> > fListOfClustersZeroFromPlane;

  uint32_t *buf;
  int fSize;
  uint32_t *cut_buf;
  int cut_fSize;
  Int_t NCH;

  int fZeroSupCut;
  int fFECID;
  int fADCChannel;
  int fAPVID;
  int fAPVKey;
  TString fIsHitMaxOrTotalADCs;
  PRDMapping * fMapping;
  
  //process data with 5 sigma cut
  GEMPedestal *ped;

  int nTimeBin;
  vector<Float_t> fPedestalNoises;

  map<Int_t, GEMHit*> fListOfHitsZero;
  map<TString, list<GEMHit*> > fListOfHitsZeroFromPlane;

};


#endif
