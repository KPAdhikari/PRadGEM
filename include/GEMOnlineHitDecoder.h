#ifndef _GEMONLINEHITDECODER_H__
#define _GEMONLINEHITDECODER_H__

#include "GEMDataStruct.h"
#include "TString.h"
#include "TH1F.h"
#include <unordered_map>
#include <map>
#include <list>
#include <vector>
#include <string>

class GEMRawDecoder;
class GEMMapping;
class GEMHit;
class GEMPedestal;
class GEMCluster;

class GEMOnlineHitDecoder
{
public:
    GEMOnlineHitDecoder();
    ~GEMOnlineHitDecoder();

    void ProcessEvent(std::unordered_map<int, std::vector<int> > &);
    void EventHandler(std::unordered_map<int, std::vector<int> > &);

    std::map< TString, std::list<GEMHit*> > GetListOfHitsFromPlanes();
    TH1F* GetHit(TString plane);

    std::map< TString, std::list<GEMHit*> > GetListOfHitsCleanFromPlanes();
    TH1F* GetCleanHit(TString plane);

    std::map<TString, std::list<GEMCluster*> > * GetListOfClustersCleanFromPlane();

    Bool_t IsADCchannelActive() ;
    void APVFindStart();
    void APVCommonModeCorrection();
    void APVCommonModeCorrectionSplit();
    void APVZeroSuppression();
    void APVEventDecoder();
    void APVEventSplitChannelsDecoder();

    // cluster
    void ComputeClusters();
    void DeleteClustersInPlaneMap() ;
    std::map < TString, std::list <GEMCluster * > > GetListOfClustersFromPlanes() { 
        return  fListOfClustersCleanFromPlane;  
    }

    TH1F* GetCluster(TString str);
    void Clear();
    void ClearClusters();
    void ClearHits();
    void SetPedestal(GEMPedestal *);
    void FillPedestal();

private:
    Int_t NCH;

    int fFECID;
    int fADCChannel;
    int fAPVID;
    int fAPVHeaderLevel;
    int fAPVKey;
    TString fIsHitMaxOrTotalADCs;
    std::vector<int> fActiveADCchannels;
    GEMMapping * mapping;
    GEMPedestal * ped;
    TString fAPVStatus;

    std::map<Int_t, GEMHit*> fListOfHits;
    std::map<Int_t, GEMHit*> fListOfHitsClean;
    std::map<TString, std::list<GEMHit*> > fListOfHitsFromPlane;
    std::map<TString, std::list<GEMHit*> > fListOfHitsCleanFromPlane;

    std::string pedestal_file;
    int fZeroSupCut;
    int fStartData;
    int nbTimeBin;
    int specialAPV;

    std::vector<int> fRawData16bits;
    std::vector <Float_t> fPedestalNoises, fPedestalOffsets;

    TString fIsClusterMaxOrTotalADCs;
    Bool_t fIsGoodClusterEvent;
    Int_t fMinClusterSize, fMaxClusterSize;

    std::map<TString, std::list<GEMCluster*> > fListOfClustersCleanFromPlane;
};


#endif
