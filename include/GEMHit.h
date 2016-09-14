/********************************************************************************
 *  GEMHit                                                                      *
 *  PRD Module Class                                                            *
 *  Author: Kondo GNANVO 12/27/2015                                             *
 *          Xinzhan Bai  03/20/2016                                             *
 *******************************************************************************/
#ifndef GEMHIT_H
#define GEMHIT_H
#include <map>
#include <TObject.h>
#include <TString.h>
#include <cassert>

class GEMMapping;

class GEMHit : public TObject
{ 
public:
    GEMHit() ;
    GEMHit(Int_t hitID, Int_t apvID, Int_t chNo, 
           Int_t zeroSupCut, TString isHitMaxOrTotalADCs) ;
    ~GEMHit() ;

    void TimingFindPeakTimeBin() ;
    Bool_t IsSortable() const { return kTRUE; }

    // sort hit according to the strip number
    Int_t Compare(const TObject *obj) const { return (fStripNo > ((GEMHit*)obj)->GetStripNo()) ? 1 : -1; }
    void ComputePosition() ;

    void AddTimeBinADCs(Int_t timebin, Float_t charges) ;
    void ClearTimeBinADCs() { fTimeBinADCs.clear() ;}

    Int_t GetAPVID() { return fAPVID;}
    Int_t GetAPVOrientation() { return fAPVOrientation;}
    Int_t GetAPVIndexOnPlane() { return fAPVIndexOnPlane;}
    Int_t GetNbAPVsFromPlane() { return fNbOfAPVsOnPlane;}
    Float_t GetHitADCs() { return fHitADCs;}
    void SetHitADCs(Float_t adc){fHitADCs = adc;}

    std::map <Int_t, Float_t> GetTimeBinADCs() { return fTimeBinADCs ; }

    Int_t StripMapping(Int_t chNo) ;
    Int_t APVchannelCorrection(Int_t chNo) ;
    Int_t PRadStripMapping(Int_t chNo) ;

    Int_t GetSignalPeakTimeBin() 
    {
	TimingFindPeakTimeBin() ;
	return fSignalPeakTimeBin;
    }

    TString GetPlane() { return fPlane; }
    Float_t GetPlaneSize() { return fPlaneSize; }

    TString GetDetector() { return fDetector; }
    TString GetDetectorType() { return fDetectorType; }
    TString GetReadoutBoard() { return fReadoutBoard; }
    TString GetHitMaxOrTotalADCs() { return fIsHitMaxOrTotalADCs ; } 

    void  SetStripNo() ;

    Int_t GetStripNo(){ return fStripNo; }
    Int_t GetAbsoluteStripNo() { return fAbsoluteStripNo; }

    Float_t GetStripPosition() 
    {
	ComputePosition() ;
	return fStripPosition;
    }

private:
    std::map <Int_t, Float_t> fTimeBinADCs;
    GEMMapping * mapping;

    Int_t fDetectorID, fPlaneID, fAPVID, fHitID, fAPVChNo, fStripNo, fAbsoluteStripNo;
    Int_t fAPVIndexOnPlane, fNbOfAPVsOnPlane, fAPVOrientation;
    Int_t fSignalPeakTimeBin, fZeroSupCut;   

    Float_t fHitADCs, fPeakADCs, fIntegratedADCs, fStripPosition, fPlaneSize,  fHitPedestalNoise;
    TString fPlane, fReadoutBoard, fDetectorType, fDetector, fIsHitMaxOrTotalADCs;

    Int_t NCH;
};

#endif
