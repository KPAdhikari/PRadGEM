/*******************************************************************************
 *  GEMCluster                                                                  *
 *  GEM Module Class                                                            *
 *  Author: Kondo GNANVO 18/08/2010                                             *
 *******************************************************************************/
#ifndef GEMCLUSTER_H
#define GEMCLUSTER_H

#include <TObject.h>
#include <TObjArray.h>

class GEMHit;

class GEMCluster: public TObject
{
public:
    GEMCluster(int minClusterSize, 
	    int maxClusterSize, 
	    TString isMaximumOrTotalCharges);

    ~GEMCluster();
    Bool_t IsSortable() const { return kTRUE; }
    TObjArray* GetArrayOfHits() { return fArrayOfHits; }

    GEMHit * GetHit(int i) 
    { 
	TObjArray &temp = * fArrayOfHits ;
	return (GEMHit *)temp[i];
    }

    void SetMinClusterSize(int min) { fMinClusterSize = min ; }
    void SetMaxClusterSize(int max) { fMaxClusterSize = max ; }
    void AddHit(GEMHit * h);
    int Compare(const TObject *obj) const ; 
    int & GetNbOfHits(){ return fNbOfHits ; }
    TString GetPlane(){	return fPlane; }
    void SetPlane(TString planename){ fPlane = planename; }
    int GetNbAPVsFromPlane(){ return fNbAPVsOnPlane; }
    void SetNbAPVsFromPlane(int nb){ fNbAPVsOnPlane = nb; }
    int GetAPVIndexOnPlane() { return fapvIndexOnPlane; }
    void SetAPVIndexOnPlane(int nb) { fapvIndexOnPlane = nb; }
    float GetPlaneSize() { return fPlaneSize; }
    void SetPlaneSize(float planesize) { fPlaneSize = planesize; }
    float & GetClusterPosition(){ return fposition; }
    float & GetClusterCentralStrip() { return fclusterCentralStrip; }

    void Timing() ;
    int GetClusterTimeBin()  ;
    int GetClusterPeakTimeBin() { return fClusterPeakTimeBin ; }
    float GetClusterADCs() ;
    void SetClusterADCs(float adc) { fClusterSumADCs = adc; }

    void ClearArrayOfHits();
    Bool_t IsGoodCluster() ;
    void ClusterCentralStrip();
    void ClusterPositionPulseHeghtWeight() ;
    std::vector< float > GetClusterTimeBinADCs () { return fClusterTimeBinADCs;  }
    void ComputeClusterPosition() ;

private:
    int fNbOfHits;
    TObjArray * fArrayOfHits;  // APV hits table
    int fClusterPeakTimeBin, fClusterTimeBin;
    float fClusterPeakADCs, fClusterTimeBinADC, fClusterSumADCs, fposition, fclusterCentralStrip, fstrip, fPlaneSize;
    int fapvID, fStripNo, fAbsoluteStripNo, fapvIndexOnPlane, fNbAPVsOnPlane, fMinClusterSize, fMaxClusterSize;   
    TString fIsClusterMaxOrSumADCs, fPlane;
    Bool_t fIsGoodCluster ;

    std::vector< float > fClusterTimeBinADCs ;
}; 
#endif
