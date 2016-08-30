//=======================================================================//
//                                                                       //
// Xinzhan Bai                                                           //
// 08/02/2016                                                            //
//=======================================================================//
#include "GEMOnlineHitDecoder.h"
#include "GEMRawDecoder.h"
#include "GEMMapping.h"
#include "GEMHit.h"
#include "GEMPedestal.h"
#include "GEMCluster.h"
#include "TCanvas.h" // Test

#define Time_Sample_SKIP 12

using namespace std;

static string trim(const string &str, const string &w = " \t\n\r")
{

    const auto strBegin = str.find_first_not_of(w);
    if (strBegin == string::npos)
	return ""; // no content

    const auto strEnd = str.find_last_not_of(w);

    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


GEMOnlineHitDecoder::GEMOnlineHitDecoder()
{
    cout<<"GEMOnlineHitDecoder Constructor..."<<endl;
    NCH = 128;
    fZeroSupCut = 5;
    fIsHitMaxOrTotalADCs = "signalPeak"; 
    //fIsHitMaxOrTotalADCs = "integratedADCs"; 
    fStartData = 0;
    nbTimeBin = 3;
    mapping = GEMMapping::GetInstance();

    // cluster
    fMinClusterSize = 1;
    fMaxClusterSize = 20;
    fIsClusterMaxOrTotalADCs = "maximumADCs";
    //fIsClusterMaxOrTotalADCs = "totalADCs";
    fIsGoodClusterEvent = kFALSE;
}

GEMOnlineHitDecoder::~GEMOnlineHitDecoder()
{
    ClearHits();
    ClearClusters();
    Clear();
    DeleteClustersInPlaneMap();
}

void GEMOnlineHitDecoder::Clear()
{
    fListOfClustersCleanFromPlane.clear();
    fListOfHits.clear();
    fListOfHitsFromPlane.clear();
    fListOfHitsClean.clear();
    fListOfHitsCleanFromPlane.clear();

    fActiveADCchannels.clear();
}

void GEMOnlineHitDecoder::ClearHits()
{
    for(auto &i: fListOfHits)
    {
	delete i.second;
    }
    fListOfHits.clear();

    for(auto &j: fListOfHitsClean)
    {
	delete j.second;
    }
    fListOfHitsClean.clear();
}

void GEMOnlineHitDecoder::ClearClusters()
{
    for(auto &i: fListOfClustersCleanFromPlane)
    {
	for(auto &j: i.second)
	{
	    delete j;
	}
	i.second.clear();
    }
    fListOfClustersCleanFromPlane.clear();
}

void GEMOnlineHitDecoder::SetPedestal(GEMPedestal * fped)
{
    ped = fped;
}

void GEMOnlineHitDecoder::ProcessEvent( unordered_map<int, vector<int> > &event)
{
    ClearHits();
    ClearClusters();
    Clear();
    DeleteClustersInPlaneMap();
    EventHandler( event);
}

void GEMOnlineHitDecoder::FillPedestal()
{
    fPedestalNoises.clear();
    fPedestalOffsets.clear();
    fPedestalNoises = ped->GetAPVNoises(fAPVID);
    fPedestalOffsets = ped->GetAPVOffsets(fAPVID);
    fAPVStatus = mapping->GetAPVstatus(fAPVID);
    string apv_status = fAPVStatus.Data();
    apv_status = trim(apv_status);
    specialAPV = 0;
    if( apv_status != "normal" )
    {
	specialAPV = 1;
    }
}

void GEMOnlineHitDecoder::EventHandler(unordered_map<int, vector<int> > & srsSingleEvent)
{
    for(auto &i: srsSingleEvent)
    {
	fAPVID = i.first;
	fADCChannel = fAPVID & 0xf;
	fFECID = (fAPVID >> 4)&0xf;
	if(IsADCchannelActive())
	{
	    fRawData16bits.clear();
	    fRawData16bits = i.second;

	    FillPedestal();

	    assert( fRawData16bits.size() > 0 );

	    if( ! specialAPV )
		APVEventDecoder();
	    else
		APVEventSplitChannelsDecoder();
	}
    }
    GetListOfHitsFromPlanes();
    GetListOfHitsCleanFromPlanes();

    // cluster
    ComputeClusters();
}

void GEMOnlineHitDecoder::APVFindStart()
{
    int idata = 0;
    int Size = fRawData16bits.size();
    for(idata = 0; idata != Size; idata++) 
    {
	if (fRawData16bits[idata] < fAPVHeaderLevel) 
	{
	    idata++ ;
	    if (fRawData16bits[idata] < fAPVHeaderLevel) 
	    {
		idata++ ;
		if (fRawData16bits[idata] < fAPVHeaderLevel) 
		{
		    idata += 10;
		    fStartData = idata ;
		    idata = Size-1 ;
		}
	    }
	}
    }
}

void GEMOnlineHitDecoder::APVCommonModeCorrection()
{
    int first[nbTimeBin];
    int last[nbTimeBin];
    first[0] = fStartData;
    last[0] = first[0] + NCH;
    for(int i=1;i<nbTimeBin;i++)
    {
	first[i] = last[i-1]+Time_Sample_SKIP;
	last[i] = first[i] +NCH;
    }
    float threshold = 20.;

    if( last[nbTimeBin-1] > fRawData16bits.size() )
    {
	cout<<"###WARNING:  GEMHitDecoder: APV Raw Data out of Range"<<endl;
	cout<<"      FEC: "<<fFECID<<" , "<<"Channel: "<<fADCChannel<<endl;
	cout<<"      please adjust your header level in mapping file"<<endl;
	cout<<"      increase recommended..."<<endl;
	cout<<"      normal APV."<<endl;
    }

    for(Int_t timebin = 0; timebin < nbTimeBin; timebin++) 
    {
	// offset correction
	for(int i=0;i<NCH;i++)
	{
	    fRawData16bits[first[timebin]+i] -= fPedestalOffsets[i];
	}

	// common mode calculation
	int nbHits = 0;
	float commonMode = 0.;
	for(int i=0; i<NCH; i++)
	{
	    if( ( fRawData16bits[first[timebin] + i] + fPedestalNoises[i]* threshold ) > 0 )
	    {
		commonMode += fRawData16bits[first[timebin] + i];
		nbHits++;
	    }
	}
	commonMode = commonMode/(float)nbHits;

	// common mode correction
	for(int i=0;i<NCH;i++)
	{
	    fRawData16bits[first[timebin]+i]-=commonMode;
	}
    }
}

void GEMOnlineHitDecoder::APVZeroSuppression()
{
    int first[nbTimeBin];
    int last[nbTimeBin];
    first[0] = fStartData;
    last[0] = first[0] + NCH;
    for(int i=1;i<nbTimeBin;i++)
    {
	first[i] = last[i-1]+Time_Sample_SKIP;
	last[i] = first[i] +NCH;
    }

    for(int i=0;i<NCH;i++)
    {
	float average = 0.;
	for(int j=0;j<nbTimeBin;j++)
	{
	    average += -fRawData16bits[ first[j] +i];
	}
	average/=(float)nbTimeBin;

	if(fZeroSupCut != 0)
	{
	    if(average > (fZeroSupCut * fPedestalNoises[i]))
	    {
		Int_t hitID = (fAPVKey << 8) | i;
		if(!fListOfHitsClean[hitID])
		{
		    GEMHit * hit = new GEMHit(hitID, fAPVID, i, fZeroSupCut, fIsHitMaxOrTotalADCs) ;
		    fListOfHitsClean[hitID] = hit ;
		}
		for(int k=0;k<nbTimeBin;k++)
		{
		    float data = - fRawData16bits[ first[k] + i];; 
		    fListOfHitsClean[hitID]->AddTimeBinADCs(k, data) ;
		}
	    }
	}
	else
	{
	    Int_t hitID = (fAPVKey << 8) | i;
	    if(!fListOfHits[hitID])
	    {
		GEMHit * hit = new GEMHit(hitID, fAPVID, i, fZeroSupCut, fIsHitMaxOrTotalADCs) ;
		fListOfHits[hitID] = hit ;
	    }
	    for(int k=0;k<nbTimeBin;k++)
	    {
		float data = - fRawData16bits[ first[k] + i];; 
		fListOfHits[hitID]->AddTimeBinADCs(k, data) ;
	    }
	}
    }
}

void GEMOnlineHitDecoder::APVCommonModeCorrectionSplit()
{
    int first[nbTimeBin];
    int last[nbTimeBin];
    first[0] = fStartData;
    last[0] = first[0] + NCH;
    for(int i=1;i<nbTimeBin;i++)
    {
	first[i] = last[i-1]+Time_Sample_SKIP;
	last[i] = first[i] +NCH;
    }
    int thresh1 = 200;
    int thresh2 = 20;

    if(last[nbTimeBin-1] > fRawData16bits.size())
    {
	cout<<"###WARNING:  GEMHitDecoder: APV Raw Data out of Range"<<endl;
	cout<<"      FEC: "<<fFECID<<"  , "<<"  Channel: "<<fADCChannel<<endl;
	cout<<"      please adjust your header level in mapping file"<<endl;
	cout<<"      increase recommended..."<<endl;
	cout<<"      special APV."<<endl;
    }

    for(int i=0;i<nbTimeBin;i++)
    {
	//offset correction
	for(int chNo=0;chNo<NCH;chNo++)
	{
	    fRawData16bits[first[i]+chNo] -= fPedestalOffsets[chNo];
	}
	//common mode calculation
	float commonMode1 = 0.;
	float commonMode2 = 0.;
	int n1 = 0;
	int n2 = 0;
	for(int chNo=0;chNo<NCH;chNo++)
	{
	    if(mapping->GetPRadStripMapping(fAPVID, chNo)<16)
	    {
		if( (fRawData16bits[first[i] + chNo] + thresh1 * fPedestalNoises[chNo]) > 0 )
		{
		    commonMode1 += fRawData16bits[first[i] + chNo];
		    n1++;
		}
	    }
	    else
	    {
		if( (fRawData16bits[first[i] + chNo] + thresh2 * fPedestalNoises[chNo]) > 0 )
		{
		    commonMode2 += fRawData16bits[first[i] + chNo];
		    n2++;
		}
	    }
	}
	commonMode1 = commonMode1/(float)n1;
	commonMode2 = commonMode2/(float)n2;

	//common mode correction
	for(int chNo = 0; chNo<NCH;chNo++)
	{
	    if(mapping->GetPRadStripMapping(fAPVID, chNo)<16)
	    {
		fRawData16bits[ first[i] + chNo] -= commonMode1;
	    }
	    else
	    {
		fRawData16bits[ first[i] + chNo] -= commonMode2;
	    }
	}
    }
}

void GEMOnlineHitDecoder::APVEventDecoder() 
{
    APVFindStart();
    APVCommonModeCorrection();
    APVZeroSuppression();
}

void GEMOnlineHitDecoder::APVEventSplitChannelsDecoder() 
{
    APVFindStart();
    APVCommonModeCorrectionSplit();
    APVZeroSuppression();
}

map < TString, list <GEMHit * > > GEMOnlineHitDecoder::GetListOfHitsFromPlanes() 
{
    map < Int_t, GEMHit * >::const_iterator hit_itr ;
    for(hit_itr = fListOfHits.begin(); hit_itr != fListOfHits.end(); ++hit_itr) 
    { 
	GEMHit * hit = (* hit_itr).second ;
	TString planename = hit->GetPlane() ;
	fListOfHitsFromPlane[planename].push_back(hit) ;
    }
    return fListOfHitsFromPlane ;
}

TH1F* GEMOnlineHitDecoder::GetHit(TString str)
{
    TH1F * h1;
    int NN = mapping->GetNbOfPlane();
    int nbDetector = mapping->GetNbOfDetectors();

    for(int i=0;i<nbDetector;i++)
    {
	TString detectorName = mapping->GetDetectorFromID(i);
	list<TString> planeList = mapping->GetPlaneListFromDetector(detectorName);
	int nbPlane = planeList.size();

	list<TString>::iterator it;
	for(it=planeList.begin();it!=planeList.end();++it)
	{ 
	    if(*it == str)
	    { 
		TString hh = detectorName+"_"+(*it)+"_hit_distribution";
		int N = (int) (mapping->GetPlaneSize(*it))/0.4;
		h1 = new TH1F(hh, hh, N, -mapping->GetPlaneSize(*it)/2, mapping->GetPlaneSize(*it)/2 );
		list< GEMHit* > hitList = fListOfHitsFromPlane[ *it  ];
		list< GEMHit* >::iterator hit_it;
		for(hit_it=hitList.begin(); hit_it!=hitList.end();++hit_it)
		{
		    Float_t pos = (*hit_it) -> GetStripPosition();
		    Float_t adc = (*hit_it) -> GetHitADCs(); 
		    h1 -> Fill(pos, adc);
		}
	    }
	}
    }
    return h1;
}

map < TString, list <GEMHit * > > GEMOnlineHitDecoder::GetListOfHitsCleanFromPlanes() 
{
    map < Int_t, GEMHit * >::const_iterator hit_itr ;
    for(hit_itr = fListOfHitsClean.begin(); hit_itr != fListOfHitsClean.end(); ++hit_itr) 
    { 
	GEMHit * hit = (* hit_itr).second ;
	TString planename = hit->GetPlane() ;
	fListOfHitsCleanFromPlane[planename].push_back(hit) ;
    }
    return fListOfHitsCleanFromPlane ;
}

TH1F* GEMOnlineHitDecoder::GetCleanHit(TString str)
{
    TH1F * h2;

    int NN = mapping->GetNbOfPlane();
    int nbDetector = mapping->GetNbOfDetectors();

    for(int i=0;i<nbDetector;i++)
    {
	TString detectorName = mapping->GetDetectorFromID(i);
	list<TString> planeList = mapping->GetPlaneListFromDetector(detectorName);
	int nbPlane = planeList.size();

	list<TString>::iterator it;
	for(it=planeList.begin();it!=planeList.end();++it)
	{
	    if(*it == str)
	    {
		TString hh = detectorName+"_clean_"+(*it)+"_hit_distribution";
		h2 = new TH1F(hh, hh, 10000, -mapping->GetPlaneSize(*it)/2-100, mapping->GetPlaneSize(*it)/2+100 );
		list< GEMHit* > hitList = fListOfHitsCleanFromPlane[ *it  ];
		list< GEMHit* >::iterator hit_it;
		for(hit_it=hitList.begin(); hit_it!=hitList.end();++hit_it)
		{
		    Float_t pos = (*hit_it) -> GetStripPosition();
		    Float_t adc = (*hit_it) -> GetHitADCs(); 
		    h2 -> Fill(pos, adc);
		}
	    }
	}
    }
    return h2;
}

Bool_t GEMOnlineHitDecoder::IsADCchannelActive() 
{
    Bool_t isADCchannelActive = kFALSE ;
    fActiveADCchannels = mapping->GetActiveADCchannels(fFECID) ;
    if ( find ( fActiveADCchannels.begin(), fActiveADCchannels.end(), fADCChannel ) != fActiveADCchannels.end() ) 
    {
	fAPVKey            = mapping->GetAPVNoFromID(fAPVID);
	fAPVHeaderLevel    = mapping->GetAPVHeaderLevelFromID(fAPVID);
	isADCchannelActive = kTRUE ;
    }
    return isADCchannelActive ;
}

//===========================================================================//
//                            Clusters Information                           //
//===========================================================================//
static Bool_t CompareStripNo( TObject *obj1, TObject *obj2) {
    Bool_t compare ;
    if ( ( (GEMHit*) obj1 )->GetStripNo() < ( ( GEMHit*) obj2 )->GetStripNo() ) compare = kTRUE ;
    else compare = kFALSE ;
    return compare ;
}

static Bool_t CompareHitADCs( TObject *obj1, TObject *obj2) {
    Bool_t compare ;
    if ( ( (GEMHit*) obj1 )->GetHitADCs() > ( ( GEMHit*) obj2 )->GetHitADCs()) compare = kTRUE ;
    else compare = kFALSE ;
    return compare ;
}

static Bool_t CompareClusterADCs( TObject *obj1, TObject *obj2) {
    Bool_t compare ;
    if ( ( (GEMCluster*) obj1 )->GetClusterADCs() > ( ( GEMCluster*) obj2 )->GetClusterADCs()) compare = kTRUE ;
    else compare = kFALSE ;
    return compare ;
}

#define NOISE_SIGMA 14
void GEMOnlineHitDecoder::ComputeClusters()
{
    map< TString, list<GEMHit*> >::iterator plane_hits_it;

    for(plane_hits_it=fListOfHitsCleanFromPlane.begin(); plane_hits_it!=fListOfHitsCleanFromPlane.end(); ++plane_hits_it)
    {
	TString plane = (*plane_hits_it).first;
	list<GEMHit*> hitsFromPlane = (*plane_hits_it).second;
	hitsFromPlane.sort(CompareStripNo);
	Int_t listSize = hitsFromPlane.size();

	if( listSize < fMinClusterSize){
	    fIsGoodClusterEvent = kFALSE;
	    continue;
	}

	ComputePlaneCluster(plane, hitsFromPlane);

	fListOfClustersCleanFromPlane[plane].sort(CompareClusterADCs);
	hitsFromPlane.clear();
    }
}

void GEMOnlineHitDecoder::ComputePlaneCluster(TString &plane, list<GEMHit*> &hitsFromPlane)
{
    list<GEMHit*>::iterator hit_it = hitsFromPlane.begin();
    list<GEMHit*>::iterator next = hit_it;

    for(hit_it;hit_it!=hitsFromPlane.end();++hit_it)
    {
	++next;
	// remove first 16 strips( apv index 0 on X side) and last 16 strips (apv index 10 on X side)
	if( plane.Contains("X")){
	    if( (*hit_it)->GetStripNo() < 16 )
		continue;
	    if( (*hit_it)->GetStripNo() > 1391 )
		continue;
	    if( (*next)->GetStripNo()>1391)
		next = hitsFromPlane.end();
	}

	Bool_t extremum = kFALSE;

	GEMCluster *new_cluster = new GEMCluster(fMinClusterSize, fMaxClusterSize, fIsClusterMaxOrTotalADCs);
	new_cluster->SetNbAPVsFromPlane((*hit_it)->GetNbAPVsFromPlane());
	new_cluster->SetAPVIndexOnPlane((*hit_it)->GetAPVIndexOnPlane());
	new_cluster->SetPlaneSize((*hit_it)->GetPlaneSize());
	new_cluster->SetPlane((*hit_it)->GetPlane());
	new_cluster->AddHit( (*hit_it) );

	// for single hit cluster
	if( (*hit_it)->GetStripNo() - (*next)->GetStripNo() == -1 ){
	    ++hit_it;
	    ++next;
	}

	if( hit_it != hitsFromPlane.end())
	{
	    while( (*hit_it)->GetStripNo() - (*next)->GetStripNo() == -1)
	    {
		if(extremum){
		    if( (*next)->GetHitADCs()  - (*hit_it)->GetHitADCs() > NOISE_SIGMA){
			(*hit_it)->SetHitADCs( (*hit_it)->GetHitADCs()/2);
			new_cluster->AddHit( (*hit_it) );
			--hit_it;
			--next;
			break;
		    }
		}
		else {
		    if( (*hit_it)->GetHitADCs() - (*next)->GetHitADCs() > NOISE_SIGMA)
			extremum = kTRUE;
		}

		new_cluster -> AddHit( (*hit_it) );
		++hit_it;
		++next;

		if( hit_it == hitsFromPlane.end() )
		    break;
	    }

	    if( (*hit_it)->GetStripNo() - (*next)->GetStripNo() != -1)
		new_cluster -> AddHit(*hit_it);
	}
	if(new_cluster->IsGoodCluster()){
	    new_cluster->ComputeClusterPosition();
	    fListOfClustersCleanFromPlane[plane].push_back(new_cluster);
	}
	else
	    delete new_cluster;
    }
}

TH1F* GEMOnlineHitDecoder::GetCluster(TString str)
{
    TH1F * hc1;

    int NN = mapping->GetNbOfPlane();

    int nbDetector = mapping->GetNbOfDetectors();

    for(int i=0;i<nbDetector;i++)
    {
	TString detectorName = mapping->GetDetectorFromID(i);
	list<TString> planeList = mapping->GetPlaneListFromDetector(detectorName);
	int nbPlane = planeList.size();

	list<TString>::iterator it;
	for(it=planeList.begin();it!=planeList.end();++it)
	{
	    if(*it == str)
	    {
		TString hh = detectorName+"_"+(*it)+"_Cluster_Distribution";
		hc1 = new TH1F(hh, hh, 10000, -mapping->GetPlaneSize(*it)/2-100, mapping->GetPlaneSize(*it)/2+100 );
		list< GEMCluster* > clusterList = fListOfClustersCleanFromPlane[ *it  ];
		list< GEMCluster* >::iterator cluster_it;
		for(cluster_it=clusterList.begin(); cluster_it!=clusterList.end();++cluster_it)
		{
		    Float_t pos = (*cluster_it) -> GetClusterPosition();
		    Float_t adc = (*cluster_it) -> GetClusterADCs(); 
		    hc1 -> Fill(pos, adc);
		}
	    }
	}
    }
    return hc1;
}

void GEMOnlineHitDecoder::DeleteClustersInPlaneMap() 
{
    map < TString, list <GEMCluster *> >::const_iterator itr ;
    for (itr = fListOfClustersCleanFromPlane.begin(); itr != fListOfClustersCleanFromPlane.end(); itr++) 
    {
	list <GEMCluster *> listOfCluster = (*itr).second ;
	listOfCluster.clear() ;
    }
    fListOfClustersCleanFromPlane.clear() ;
}

map<TString, std::list<GEMCluster*> > * GEMOnlineHitDecoder::GetListOfClustersCleanFromPlane()
{
    return &fListOfClustersCleanFromPlane;
}
