#include "GEMHit.h"
//ClassImp(GEMHit);

GEMHit::GEMHit(Int_t hitID, Int_t apvID, Int_t chNo, Int_t zeroSupCut, TString isHitMaxOrTotalADCs) {
  
  fTimeBinADCs.clear();
  fHitADCs = -10000 ;
  fPeakADCs = 0.0 ;
  fIntegratedADCs = 0.0 ;
  fSignalPeakTimeBin = 0 ;
 
  NCH = 128;

  fIsHitMaxOrTotalADCs = isHitMaxOrTotalADCs ;
  fZeroSupCut = zeroSupCut ;

  PRDMapping * mapping = PRDMapping::GetInstance() ;
  fAPVID             = apvID ;
  fAPVIndexOnPlane   = mapping->GetAPVIndexOnPlane(fAPVID);
  fAPVOrientation    = mapping->GetAPVOrientation(fAPVID);
  fPlane             = mapping->GetPlaneFromAPVID(fAPVID); 
  fPlaneID           = mapping->GetPlaneID(fPlane) ;
  fDetector          = mapping->GetDetectorFromPlane(fPlane) ;
  fDetectorID        = mapping->GetDetectorID(fDetector) ;
  fDetectorType      = mapping->GetDetectorTypeFromDetector(fDetector) ;
  fReadoutBoard      = mapping->GetReadoutBoardFromDetector(fDetector) ;
  fPlaneSize         = mapping->GetPlaneSize(fPlane);
  fNbOfAPVsOnPlane   = mapping->GetNbOfAPVsOnPlane(fPlane);
  fAPVChNo = chNo ;
  SetStripNo() ;
  ComputePosition() ;

}

GEMHit::~GEMHit(){
  fTimeBinADCs.clear() ;
}

void GEMHit::TimingFindPeakTimeBin() {
  Bool_t timingStatus = kTRUE ;
  Float_t currentMax = 0.0 ;
  map <Int_t, Float_t>::const_iterator  max_itr ;
  for(max_itr = fTimeBinADCs.begin(); max_itr != fTimeBinADCs.end(); ++max_itr) { 
    if (max_itr->second > currentMax) {
      currentMax = max_itr->second ;
      fSignalPeakTimeBin = max_itr->first;
    }
  }
}

void GEMHit::AddTimeBinADCs(Int_t timebin, Float_t charges) {
  assert(timebin>=0);
  fTimeBinADCs[timebin] = charges;
  if(fZeroSupCut > 0 ) {
    //if (charges >  fHitADCs) fPeakADCs = charges ;
    if (charges >  fPeakADCs) fPeakADCs = charges ; // fix xb
    fIntegratedADCs += charges ;
    if (fIsHitMaxOrTotalADCs == "integratedADCs") fHitADCs = fIntegratedADCs ;
    else                                          fHitADCs =  fPeakADCs;
  }

  else {
    if (timebin ==0)
      fHitADCs = charges ; 
    else 
    {
      fHitADCs += charges;
      fHitADCs/=(timebin+1);
    }
  }
}

void GEMHit::ComputePosition() {

  Float_t pitch = 0.4; //0.4mm
  if ( (fDetectorType == "PRADGEM") && (fPlane.Contains("X")) ) {
    
    //mapping fix by xb
    fStripPosition =  -0.5 * (fPlaneSize + 12.8 /* 32 strips*0.4 pitch */-pitch ) + (pitch * fStripNo) ;

  }
  else {
    fStripPosition =  -0.5 * (fPlaneSize - pitch) + (pitch * fStripNo) ;
  }
}

void GEMHit::SetStripNo() {

  fAbsoluteStripNo = StripMapping(fAPVChNo) ;

  if ( (fDetectorType == "PRADGEM")  && (fPlane.Contains("X")) ) 
    {
      Int_t nbAPVsOnPlane = fNbOfAPVsOnPlane - 1 ;
      if(fAPVIndexOnPlane != 11)
      {
        if (fAPVIndexOnPlane > fNbOfAPVsOnPlane) fStripNo = -100000000 ;
        if(fAPVOrientation == 0) fAbsoluteStripNo = 127 - fAbsoluteStripNo ;
        //fStripNo = (fAbsoluteStripNo-16) + (NCH * (fAPVIndexOnPlane % (nbAPVsOnPlane))) ; //this will cause chamber center to shift, do not use
        fStripNo = (fAbsoluteStripNo) + (NCH * (fAPVIndexOnPlane % (nbAPVsOnPlane))) ; // take chamber center as origin, should not shift
      }
      else if (fAPVIndexOnPlane == 11) 
	{
          
	  Int_t apvIndexOnPlane = fAPVIndexOnPlane - 1 ;
	  if (apvIndexOnPlane > fNbOfAPVsOnPlane) fStripNo = -100000000 ;
	  if(fAPVOrientation == 0) fAbsoluteStripNo = 127 - fAbsoluteStripNo ;
	  //fStripNo = (fAbsoluteStripNo - 32 ) + (NCH * (apvIndexOnPlane % nbAPVsOnPlane)) ;
	  fStripNo = (fAbsoluteStripNo - 16) + (NCH * (apvIndexOnPlane % nbAPVsOnPlane)) ;
          
	}
    }

  else 
    {
  
      if (fAPVIndexOnPlane > fNbOfAPVsOnPlane) fStripNo = -100000000 ;
      if(fAPVOrientation == 0) fAbsoluteStripNo = 127 - fAbsoluteStripNo ;
      fStripNo = fAbsoluteStripNo + (NCH * (fAPVIndexOnPlane % fNbOfAPVsOnPlane)) ;
  
    }
}

Int_t GEMHit::StripMapping(Int_t chNo) {
  chNo = APVchannelCorrection(chNo) ;
  if (fDetectorType == "PRADGEM") chNo = PRadStripMapping(chNo) ;
  return chNo ;
}

Int_t GEMHit::APVchannelCorrection(Int_t chNo) {
  chNo = (32 * (chNo%4)) + (8 * (Int_t)(chNo/4)) - (31 * (Int_t)(chNo/16)) ;
  return chNo ;
}

Int_t GEMHit::PRadStripMapping(Int_t chNo) {
  if ( (fDetectorType == "PRADGEM") && (fPlane.Contains("X")) && (fAPVIndexOnPlane == 11) ) {
    if (chNo % 2 == 0) chNo = ( chNo / 2) + 48 ; // originally was 48 xb: very important
    else
      if (chNo < 96) chNo = (95 - chNo) / 2 ;
      else           chNo = 127 + (97 - chNo) / 2 ;
  }
  else { // NON (fDetectorType == "PRADGEM") && (fPlane.Contains("Y")) && (fAPVIndexOnPlane == 11)
    if (chNo % 2 == 0) chNo = ( chNo / 2) + 32 ;
    else 
      if (chNo < 64) chNo = (63 - chNo) / 2 ;
      else           chNo = 127 + (65 - chNo) / 2 ;
  }
  return chNo ;
}
