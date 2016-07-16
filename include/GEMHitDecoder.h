#ifndef __GEMHITDECODER__
#define __GEMHITDECODER__
/*******************************************************************************
 *  AMORE FOR PRD - PRD                                                         *
 *  GEMHitDecoder                                                               *
 *  PRD Module Class                                                            *
 *  Author: Kondo GNANVO 12/27/2015                                             *
 *          Xinzhan Bai  03/22/2016                                             *
 *******************************************************************************/

//#define NCH 128  
//#define PI 3.14159265359 
//#define CEILING 4096

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <sstream>
#include <fstream>
#include <iostream>
#include "TObject.h"
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TColor.h>
#include <TCanvas.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TList.h"
#include "TString.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMath.h"

#include "TTree.h"

#include <map>
#include <list>
#include <vector>
#include <valarray>
#include <numeric>

#include "evioUtil.hxx"
#include "evioFileChannel.hxx"

#include "GEMHit.h"
#include "GEMPedestal.h"
#include "PRDMapping.h"
//#include "GEMRawDecoder.h"
#endif

using namespace std;
using namespace evio;

//class GEMHit;

//class GEMHitDecoder : public TObject {
class GEMHitDecoder {

 public:

  //GEMHitDecoder(Int_t zeroSup); //zeroSup will take from configure file
  GEMHitDecoder(string);
  GEMHitDecoder();
  ~GEMHitDecoder(); 

  void Clear()  ;
  void ClearListOfHits()  ;
  void ClearListOfHitsFromPlane() ;

  //void FECEventDecoder(GEMPedestal * ped, UInt_t * buf) ;
  Int_t ProcessAllEvents(Int_t evt=-1);
  void EventHandler(map<int, map<int, vector<int> > > , GEMPedestal*);

  void APVEventDecoder() ;
  void APVEventSplitChannelsDecoder();
  //void ComputeHits(GEMPedestal * ped) ;
  //void Convert32bitTo16bit(UInt_t rawData16bits) ;

  map < Int_t, GEMHit * >  GetListOfHits() {return fListOfHits;}
  map < TString, list <GEMHit * > >  GetListOfHitsFromPlanes() ;

  Int_t GetEventNumber() {return fEventNb;}

  Bool_t IsADCchannelActive() ;

 private:

  PRDMapping * fMapping ;
  map<int, map<int, vector<int> > > mSrsSingleEvent;

  Int_t fNbOfChannels, fEventNb, fZeroSupCut;
  Int_t fNbOfAPVs, fFECID, fCurrentFEC, fADCChannel, fAPVID, fAPVKey, fDetectorID, fPlaneID;
  Int_t fAPVIndexOnPlane, fAPVOrientation, fNbOfAPVsOnPlane;
  Int_t fStartData;
  TString fRunName, fRunType, fAPVName, fPlane, fDetector, fDetectorType, fReadoutBoard, fIsHitMaxOrTotalADCs ;

  TString fAPVStatus;  // xb

  Float_t fPlaneSize, fEtaSectorPos, fAPVGain, fMeanAPVnoise;
  Bool_t  fIsNewFECdataFlag;
  UInt_t fAPVHeaderLevel ;

  vector <Int_t>  fActiveADCchannels ;
  //vector<UInt_t>  fRawData16bits;
  vector<int> fRawData16bits;

  vector <Float_t> fPedestalNoises, fPedestalOffsets;
  vector <Float_t> fPedestalNoises_1stSet, fPedestalNoises_2ndSet, fPedestalOffsets_1stSet, fPedestalOffsets_2ndSet;
  map < Int_t, GEMHit * > fListOfHits ;
  map < TString, list <GEMHit * > > fListOfHitsFromPlane ;

private:
  // xb 
  Int_t NCH;
  GEMPedestal *ped;

  string file_name;
  string pedestal_file;
  set<int> FECs; //FEC ID

public:
  // xb for debug
  void ShowHits();

  //ClassDef(GEMHitDecoder,3);
};

#endif
