#include "GEMCluster.h"
#include "GEMHit.h"
#include <TH1F.h>
#include <map>

GEMCluster::GEMCluster(Int_t minClusterSize, Int_t maxClusterSize, TString isMaximumOrTotalADCs) 
{
    fNbOfHits = 0;  
    fstrip = 0;
    fClusterSumADCs = 0; 
    fClusterPeakADCs = 0;
    fClusterTimeBinADC = 0 ;
    fClusterTimeBin = 0 ;
    fClusterPeakTimeBin = 0 ;
    fposition = 0;
    fclusterCentralStrip = 0;
    fMinClusterSize = minClusterSize;
    fMaxClusterSize = maxClusterSize;
    fIsClusterMaxOrSumADCs = isMaximumOrTotalADCs ;
    fPlaneSize = 512. ;
    fPlane = "GEM1X" ;
    fNbAPVsOnPlane = 10 ;
    fArrayOfHits = new TObjArray(maxClusterSize);
    fIsGoodCluster = kTRUE ;
}

GEMCluster::~GEMCluster() 
{
    fArrayOfHits->Clear();
    delete fArrayOfHits;
    fClusterTimeBinADCs.clear() ;
}

void GEMCluster::Timing() 
{
    Int_t nBins = fClusterTimeBinADCs.size() ;
    TH1F * timeBinHist = new TH1F("timeBinHist", "timeBinHist", nBins, 0, (nBins-1) ) ;
    for (Int_t k = 0; k < nBins; k++)  timeBinHist->Fill(k,fClusterTimeBinADCs[k]) ;
    fClusterTimeBin = timeBinHist->GetMaximumBin() ; 
    delete timeBinHist ;
}

Int_t GEMCluster::GetClusterTimeBin() 
{  
    float q;
    TObjArray &temp = *fArrayOfHits;

    Int_t nbofhits =  GetNbOfHits() ;

    fClusterTimeBinADC = 0 ;  // to be removed...???
    for (int i = 0; i < nbofhits; i++) {	
	std::map <Int_t, Float_t>  timeBinADCs = ((GEMHit*)temp[i])->GetTimeBinADCs() ;;
	Int_t nbOfTimeBins = timeBinADCs.size() ;
	fClusterTimeBinADCs.resize(nbOfTimeBins) ;
	for (int k = 0; k < nbOfTimeBins; k++)  fClusterTimeBinADCs[k] += timeBinADCs[k] ;
	timeBinADCs.clear() ;
    }
    Timing() ;
    return fClusterTimeBin ;
}

void GEMCluster::ClusterPositionPulseHeghtWeight() 
{  
    // Calculate the fposition and the total fClusterSumADCs
    float hitposition, q;
    TObjArray &temp = *fArrayOfHits;
    Int_t nbofhits =  GetNbOfHits() ;
    for (int i = 0; i < nbofhits; i++) {	
	q  = ((GEMHit*)temp[i])->GetHitADCs() ;
	hitposition = ((GEMHit*)temp[i])->GetStripPosition() ;
	fClusterSumADCs += q ;
	fposition += q * hitposition ;
	if (q > fClusterPeakADCs) {
	    fClusterPeakTimeBin = ((GEMHit*)temp[i])->GetSignalPeakTimeBin() ;
	    fClusterPeakADCs = q ;
	}
    }
    fposition /= fClusterSumADCs;
}

void GEMCluster::ClusterCentralStrip() 
{   
    float p, dp ;
    float dpmin = 99;
    TObjArray &temp = *fArrayOfHits;
    Int_t nbofhits =  GetNbOfHits() ;
    for (int i = 0; i < nbofhits; i++) {
	p  = ((GEMHit*)temp[i])->GetStripPosition();
	dp = fabs(fposition - p); 
	if (dp <= dpmin) {
	    fclusterCentralStrip = p;
	    dpmin = dp;
	}
    }
}

Float_t GEMCluster::GetClusterADCs() 
{
    Float_t adcs = 0 ;
    if (fIsClusterMaxOrSumADCs == "maximumADCs") {
	adcs = fClusterPeakADCs ;
    }
    else {
	adcs = fClusterSumADCs ;
    } 
    return adcs ;
}

void GEMCluster::AddHit(GEMHit *h) 
{
    fArrayOfHits->AddLast(h); 
}

void GEMCluster::ClearArrayOfHits() 
{
    fArrayOfHits->Clear();
}

Bool_t GEMCluster::IsGoodCluster() 
{
    fIsGoodCluster = kTRUE ;
    fNbOfHits = fArrayOfHits->GetEntries() ;
    if ( (fNbOfHits > fMaxClusterSize) || (fNbOfHits < fMinClusterSize) ) {
	ClearArrayOfHits() ;
	fIsGoodCluster = kFALSE ;
	fNbOfHits = fArrayOfHits->GetEntries() ;
    } 
    return fIsGoodCluster ;
}

int GEMCluster::Compare(const TObject *obj) const 
{
    int compare = (fClusterSumADCs < ((GEMCluster*)obj)->GetClusterADCs()) ? 1 : -1;
    return compare ;
}

void GEMCluster::ComputeClusterPosition() 
{
    ClusterPositionPulseHeghtWeight() ;
    ClusterCentralStrip() ;
}

