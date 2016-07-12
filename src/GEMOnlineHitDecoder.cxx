#include "GEMOnlineHitDecoder.h"

static string trim(const string &str, const string &w = " \t\n\r")
{

    const auto strBegin = str.find_first_not_of(w);
    if (strBegin == string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(w);

    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


GEMOnlineHitDecoder::GEMOnlineHitDecoder(uint32_t *rbuf, int Size, GEMPedestal *pPed)
{
  NCH = 128;
  
  fZeroSupCut = 5;
  fSaveRawHit = 0; // for online show raw hits
  fIsHitMaxOrTotalADCs = "signalPeak"; 
  fStartData = 0;
  fTimeSample = 3;
  buf = rbuf;
  fSize = Size;

  fMapping = PRDMapping::GetInstance();
  
  ped = pPed;

  FECs.clear();
  mSrsSingleEvent.clear();
  FECs = fMapping->GetBankIDSet();
 
  // cluster information  
  fMinClusterSize = 1;
  fMaxClusterSize = 20;
  //fIsClusterMaxOrTotalADCs = "totalADCs";
  fIsClusterMaxOrTotalADCs = "maximumADCs";
  fIsGoodClusterEvent = kFALSE;
  fListOfClustersCleanFromPlane.clear();

  ProcessEvent();
}

GEMOnlineHitDecoder::~GEMOnlineHitDecoder()
{
  //clear hits
  map<int, GEMHit *>::iterator it = fListOfHits.begin();
  for(;it!=fListOfHits.end();++it)
  {
    delete (it->second);
  }
  it = fListOfHitsClean.begin();
  for(;it!=fListOfHitsClean.end();++it)
  {
    delete(it->second);
  }

  //clear clusters
  if( fListOfClustersCleanFromPlane.size() > 0)
  {
    map<TString, list<GEMCluster*> >::iterator itt = fListOfClustersCleanFromPlane.begin();
    for(;itt!=fListOfClustersCleanFromPlane.end();++itt)
    {
      list<GEMCluster*>::iterator itc = (itt->second).begin();
      for(;itc!=(itt->second).end();++itc)
        delete *itc;
      (itt->second).clear();
    }
    
  }
  fListOfClustersCleanFromPlane.clear();
  fListOfHits.clear();
  fListOfHitsFromPlane.clear();
  fListOfHitsClean.clear();
  fListOfHitsCleanFromPlane.clear();
 
  fActiveADCchannels.clear();

  FECs.clear();
  mSrsSingleEvent.clear();
  DeleteClustersInPlaneMap();
}

void GEMOnlineHitDecoder::ProcessEvent()
{
  GEMRawDecoder raw_decoder(buf, fSize);
  EventHandler( raw_decoder.GetDecoded(), ped);
}

void GEMOnlineHitDecoder::EventHandler( map<int, map<int, vector<int> > > srsSingleEvent , GEMPedestal * ped )
{

  mSrsSingleEvent = srsSingleEvent;
  map<int, map<int, vector<int> > >::iterator it;
  for(it=mSrsSingleEvent.begin(); it!=mSrsSingleEvent.end(); ++it)
  {
    fFECID = it->first;
    map<int, vector<int> > fec_data = it->second;
    map<int, vector<int> >::iterator itt;
    for(itt=fec_data.begin(); itt!=fec_data.end(); ++itt)
    {
      fADCChannel = itt->first;
      fAPVID = fFECID << 4 | fADCChannel;
      if(IsADCchannelActive())
      {
        fRawData16bits.clear();
	fRawData16bits = itt->second;

	fPedestalNoises.clear();
	fPedestalOffsets.clear();
	fPedestalNoises = ped->GetAPVNoises(fAPVID);
	fPedestalOffsets = ped->GetAPVOffsets(fAPVID);
	fPedestalNoises_1stSet.clear();
	fPedestalOffsets_1stSet.clear();
	fPedestalNoises_2ndSet.clear();
	fPedestalOffsets_2ndSet.clear();

	assert( fRawData16bits.size() > 0 );
	
	// Decode APV25 DATA INTO HITS
	fAPVStatus = fMapping->GetAPVstatus(fAPVID);
	string apv_status = fAPVStatus.Data();
	apv_status = trim(apv_status);
	if( apv_status != "normal" )
	{
	  for(Int_t chNo=0; chNo<NCH; chNo++)
	  {
	    if( fMapping->GetPRadStripMapping(fAPVID, chNo) < 16)
	    {
	      fPedestalOffsets_1stSet.push_back(fPedestalOffsets[chNo]);
	      fPedestalNoises_1stSet.push_back(fPedestalNoises[chNo]);
	    }
	    else
	    {
	      fPedestalOffsets_2ndSet.push_back(fPedestalOffsets[chNo]);
	      fPedestalNoises_2ndSet.push_back(fPedestalNoises[chNo]);
	    }
	  }
          // DECODE DATA FROM SPECIAL APVs into Hits
	  APVEventSplitChannelsDecoder();
	}
	else
	{
          APVEventDecoder();
	}
      }
    }
  }

  //clear mSrsSingleEvent
  for(it = mSrsSingleEvent.begin(); it!= mSrsSingleEvent.end(); ++it)
  {
    map<int, vector<int> > fec_data = it->second;
    map<int, vector<int> >::iterator itt;
    for(itt=fec_data.begin(); itt!=fec_data.end();++itt)
    {
      itt->second.clear();
    }
    it->second.clear();
  }
  mSrsSingleEvent.clear();
  GetListOfHitsFromPlanes();
  GetListOfHitsCleanFromPlanes();

  // cluster computing      
  ComputeClusters();
}

void GEMOnlineHitDecoder::APVEventDecoder() {
  Int_t idata = 0, firstdata = 0, lastdata = 0, nbOfTimeBin = fTimeSample;
  Int_t loop_size = 0;
  Int_t threshold = 20.0;
  
  Int_t size = fRawData16bits.size() ;
  int commonMode = 0 ;

  vector<Float_t> rawDataTS, rawDataZS;
  vector<Float_t> dataTest, commonModeOffsets ;
  rawDataTS.clear();

  // COMPUTE APV25 COMMON MODE CORRECTION
  for(idata = 0; idata < size; idata++) {
    if (fRawData16bits[idata] < fAPVHeaderLevel) {
      idata++ ;
      if (fRawData16bits[idata] < fAPVHeaderLevel) {
	idata++ ;
	if (fRawData16bits[idata] < fAPVHeaderLevel) {
	  idata += 10;
	  fStartData = idata ;
	  idata = size ;
	}
      }
    }
  }
  //==========================================================================
  // BLOCK WHERE COMMON MODE CORRECTION FOR ALL APV25 TIME BINS IS COMPUTED
  //==========================================================================
  firstdata = fStartData ;
  lastdata  = firstdata  + NCH ;

  for(Int_t timebin = 0; timebin < nbOfTimeBin; timebin++) {

    if( lastdata > size ) 
      {
	cout<<"###ERROR:  GEMHitDecoder: APV Raw Data out of Range"<<endl;
	cout<<"           please adjust your header level in mapping file"<<endl;
	cout<<"           increase recommended..."<<endl;
      }
    rawDataTS.insert(rawDataTS.end(), &fRawData16bits[firstdata], &fRawData16bits[lastdata]);
    assert(rawDataTS.size() == 128 );

    // PERFORM APV25 PEDESTAL OFFSET CORRECTION  FOR A GIVEN TIME BIN
    loop_size = rawDataTS.size();
    for(int i=0;i<loop_size;i++)
      {
	rawDataTS[i]-=fPedestalOffsets[i];
      }

    dataTest.clear();
    dataTest.insert(dataTest.begin(), rawDataTS.begin(), rawDataTS.end() );
    int nbHits = 0;
    assert( rawDataTS.size() == 128 );
    for(int i=0; i<128; i++)
      {
	if( ( dataTest[i]-fPedestalNoises[i]* threshold )>=0 )
	  {
	    dataTest[i] = 0;
	    nbHits++;
	  }
      }
    assert( nbHits <= 128);

    // COMPUTE COMMON MODE FOR A GIVEN APV AND TIME BIN
    commonMode = 0.0;
    for(int i=0;i<loop_size;i++)
      {
	commonMode+=dataTest[i]; 
      }
    commonMode/=(NCH-nbHits);
   
    commonModeOffsets.push_back(commonMode) ;

    // PERFORM COMMON MODE CORRECTION FOR A GIVEN TIME BIN
    for(int i=0;i<loop_size;i++)
      {
	rawDataTS[i]-=commonMode;
      }

    //  ADC SUM OVER ALL TIME BINS USE AS THE TEST CRITERIA FOR ZERO SUPPRESSION
    if(timebin == 0)  rawDataZS.resize(rawDataTS.size()) ;
    for(int i=0;i<loop_size;i++)
      {
	rawDataZS[i]+=rawDataTS[i];
      }

    // PROCEED TO NEXT TIME BIN
    firstdata = lastdata + 12 ;
    lastdata = firstdata + NCH ;

    // CLEAR EVERYTHING
    rawDataTS.clear() ;
    dataTest.clear() ;
  }

  // ADC AVERAGE OVER ALL TIME BINS USE AS THE TEST CRITERIA FOR ZERO SUPPRESSION
  loop_size = rawDataZS.size();
  for(int i=0;i<loop_size;i++)
    {
      rawDataZS[i]/=(float)nbOfTimeBin; 
    }

  //=========================================
  // BLOCK WHERE ZERO SUPPRESSION IS COMPUTED
  //=========================================
  firstdata = fStartData ;
  lastdata  = firstdata  + NCH ;

  for(Int_t timebin = 0; timebin < nbOfTimeBin; timebin++) {

    // EXTRACT APV25 DATA FOR A GIVEN TIME BIN
    rawDataTS.insert(rawDataTS.end(), &fRawData16bits[firstdata], &fRawData16bits[lastdata]);

    for(Int_t chNo = 0; chNo < NCH; chNo++) {
      Int_t hitID = (fAPVKey << 8) | chNo ;
      Float_t data    = -(rawDataTS[chNo] - fPedestalOffsets[chNo] - commonModeOffsets[timebin]) ;      
      Float_t avgdata = -rawDataZS[chNo] ;

      // PERFORM ZERO SUPPRESSION
      if( fZeroSupCut > 0 ) 
	{
	  if( avgdata > (fZeroSupCut * fPedestalNoises[chNo]) ) 
	    {
	      if(!fListOfHitsClean[hitID])
		{
		  GEMHit * hit = new GEMHit(hitID, fAPVID, chNo, fZeroSupCut,  fIsHitMaxOrTotalADCs) ;
		  fListOfHitsClean[hitID] = hit ;
		}
	      fListOfHitsClean[hitID]->AddTimeBinADCs(timebin, data) ;
	    }
	}
     
      // NO ZERO SUPPRESSION
      if( fSaveRawHit > 0 ) 
	{
	  if(!fListOfHits[hitID]) 
	    {
	      GEMHit * hit = new GEMHit(hitID, fAPVID, chNo, 0,  fIsHitMaxOrTotalADCs) ;
	      fListOfHits[hitID] = hit ;
	    }
	  fListOfHits[hitID]->AddTimeBinADCs(timebin, data) ;
	}
      
    }    
    firstdata = lastdata + 12 ;
    lastdata = firstdata + NCH ;
    rawDataTS.clear() ;
  }
  commonModeOffsets.clear() ;
}

void GEMOnlineHitDecoder::APVEventSplitChannelsDecoder() {
  Int_t idata = 0, firstdata = 0, lastdata = 0, nbOfTimeBin = fTimeSample;
  Int_t loop_size = 0;
  Int_t threshold16 = 200.0; // track:  50 -> 200 -> 100 
  Int_t threshold112 = 20.0;
 
  Int_t size = fRawData16bits.size() ;
  Int_t size1 = fPedestalNoises_1stSet.size();
  Int_t size2 = fPedestalNoises_2ndSet.size();

  Float_t commonMode = 0 ;

  vector<Float_t> rawDataTS, rawDataTS_1stSet, rawDataTS_2ndSet; 
  vector<Float_t> rawDataZS_1stSet, rawDataZS_2ndSet;
  vector<Float_t> threshold_1stSet, threshold_2ndSet ;
  vector<Float_t> dataTest_1stSet, dataTest_2ndSet;
  vector<Float_t> commonModeOffsets_1stSet, commonModeOffsets_2ndSet ;

  rawDataTS.clear(); 
  rawDataTS_1stSet.clear();
  rawDataTS_2ndSet.clear();

  commonModeOffsets_1stSet.clear();
  commonModeOffsets_2ndSet.clear();

  // COMPUTE APV25 COMMON MODE CORRECTION
  for(idata = 0; idata < size; idata++) {
    if (fRawData16bits[idata] < fAPVHeaderLevel) {
      idata++ ;
      if (fRawData16bits[idata] < fAPVHeaderLevel) {
	idata++ ;
	if (fRawData16bits[idata] < fAPVHeaderLevel) {
	  idata += 10;
	  fStartData = idata ;
	  idata = size ;
	}
      }
    }
  }

  //=======================================================================
  // BLOCK WHERE COMMON MODE CORRECTION FOR ALL APV25 TIME BINS IS COMPUTED
  //=======================================================================
  firstdata = fStartData ;
  lastdata  = firstdata  + NCH ;

  for(Int_t timebin = 0; timebin < nbOfTimeBin; timebin++) {
    //=========================================
    // EXTRACT APV25 DATA FOR A GIVEN TIME BIN
    //=========================================
    if( lastdata > size ) 
      {
	cout<<"###Warning:  GEMHitDecoder: APV Raw Data out of Range"<<endl;
	cout<<"           please adjust your header level in mapping file"<<endl;
	cout<<"           increase recommended..."<<endl;
      }
    rawDataTS.insert(rawDataTS.end(), &fRawData16bits[firstdata], &fRawData16bits[lastdata]);
    assert( rawDataTS.size() == 128 );

    //split Channel data
    for(Int_t chNo=0; chNo<NCH;chNo++)
      {
	if(fMapping->GetPRadStripMapping(fAPVID, chNo)<16)
	  rawDataTS_1stSet.push_back(rawDataTS[chNo]);
	else
	  rawDataTS_2ndSet.push_back(rawDataTS[chNo]);
      }
    assert( rawDataTS_1stSet.size() == 16 );
    assert( rawDataTS_2ndSet.size() == 112);

    //===============================================================
    //                     1ST SET OF CHANNELS
    //===============================================================

    // PERFORM APV25 PEDESTAL OFFSET CORRECTION  FOR A GIVEN TIME BIN
    loop_size = rawDataTS_1stSet.size();
    for(int i=0;i<loop_size;i++)
      {
	rawDataTS_1stSet[i] -= fPedestalOffsets_1stSet[i];
      }

    dataTest_1stSet.clear();
    dataTest_1stSet.insert(dataTest_1stSet.begin(), rawDataTS_1stSet.begin(), rawDataTS_1stSet.end());
    assert( rawDataTS_1stSet.size() == 16 );
    int nbHits_1stSet = 0;
    for(int i=0;i<16;i++)
      {
	if( (dataTest_1stSet[i] - fPedestalNoises_1stSet[i]*threshold16) >= 0)
	  {
	    dataTest_1stSet[i] = 0;
	    nbHits_1stSet++;
	  }
      }
    assert(nbHits_1stSet < 16 );

    // COMPUTE COMMON MODE FOR A GIVEN APV AND TIME BIN
    commonMode = 0;
    for(int i=0;i<16;i++)
      {
	commonMode+=dataTest_1stSet[i];
      }
    commonMode/=(16-nbHits_1stSet);
    commonModeOffsets_1stSet.push_back(commonMode);

    // PERFORM COMMON MODE CORRECTION FOR A GIVEN TIME BIN
    loop_size = rawDataTS_1stSet.size();
    for(int i=0;i<loop_size;i++)
      {
	rawDataTS_1stSet[i]-=commonMode;
      }

    //  ADC SUM OVER ALL TIME BINS USE AS THE TEST CRITERIA FOR ZERO SUPPRESSION
    if(timebin == 0)  rawDataZS_1stSet.resize(rawDataTS_1stSet.size()) ;

    for(int i=0;i<loop_size;i++)
      {
	rawDataZS_1stSet[i] += rawDataTS_1stSet[i];
      }

    //======================================================================
    //                    2ND SET OF CHANNELS
    //======================================================================

    // PERFORM APV25 PEDESTAL OFFSET CORRECTION  FOR A GIVEN TIME BIN
    loop_size = rawDataTS_2ndSet.size();
    for(int i=0;i<loop_size;i++)
      {
	rawDataTS_2ndSet[i] -= fPedestalOffsets_2ndSet[i];
      }

    // TEST IF APV CHANNEL IS A HIT OR PEDESTAL => ROUGH TEST TO ELIMINATE HIT CHANNELS IN FOR THE COMMON MODE
    dataTest_2ndSet.clear();
    dataTest_2ndSet.insert(dataTest_2ndSet.begin(), rawDataTS_2ndSet.begin(), rawDataTS_2ndSet.end());

    int nbHits_2ndSet = 0;
    assert(size2==112);
    for(int i=0;i<size2;i++)
      {
	if( (dataTest_2ndSet[i] - fPedestalNoises_2ndSet[i]*threshold112) >= 0)
	  {
	    dataTest_2ndSet[i] = 0;
	    nbHits_2ndSet++;
	  }
      }
    assert( nbHits_2ndSet < size2);

    // COMPUTE COMMON MODE AND PERFORM COMMON MODE CORRECTION FOR A GIVEN APV AND TIME BIN FOR 2ND SET OF CHANNELS
    commonMode = 0;
    loop_size = dataTest_2ndSet.size();
    for(int i=0;i<loop_size;i++)
      {
	commonMode += dataTest_2ndSet[i];
      }
    commonMode/=(112-nbHits_2ndSet);
    commonModeOffsets_2ndSet.push_back(commonMode) ;
    for(int i=0;i<loop_size;i++)
      {
	rawDataTS_2ndSet[i] -= commonMode;
      }

    //  ADC SUM OVER ALL TIME BINS USE AS THE TEST CRITERIA FOR ZERO SUPPRESSION
    if(timebin == 0){ rawDataZS_2ndSet.clear(); rawDataZS_2ndSet.resize( rawDataTS_2ndSet.size()  );  }
    loop_size = rawDataZS_2ndSet.size();
    for(int i=0;i<loop_size;i++)
      {
	rawDataZS_2ndSet[i]+=rawDataTS_2ndSet[i];
      }
    
    //===================================
    //        CLEAR EVERYTHING
    //===================================
    rawDataTS.clear() ;
    rawDataTS_1stSet.clear() ;
    rawDataTS_2ndSet.clear() ;
    dataTest_1stSet.clear() ;
    dataTest_2ndSet.clear() ;
 
    // PROCEED TO NEXT TIME BIN
    firstdata = lastdata + 12 ;
    lastdata = firstdata + NCH ;
  }


  // ADC AVERAGE OVER ALL TIME BINS USE AS THE TEST CRITERIA FOR ZERO SUPPRESSION
  loop_size = rawDataZS_1stSet.size();
  for(int i=0;i<loop_size;i++)
    {
      rawDataZS_1stSet[i]/=(float)nbOfTimeBin;
    }
  loop_size = rawDataZS_2ndSet.size();
  for(int i=0;i<loop_size;i++)
    {
      rawDataZS_2ndSet[i]/=(float)nbOfTimeBin;
    }

  //============================================
  // BLOCK WHERE ZERO SUPPRESSION IS COMPUTED
  //============================================
  firstdata = fStartData ;
  lastdata  = firstdata  + NCH ;

  for(Int_t timebin = 0; timebin < nbOfTimeBin; timebin++) {

    // EXTRACT APV25 DATA FOR A GIVEN TIME BIN
    rawDataTS.insert(rawDataTS.end(), &fRawData16bits[firstdata], &fRawData16bits[lastdata]);

    int i1st = 0;
    int i2nd = 0;

    for(Int_t chNo = 0; chNo < NCH; chNo++) {
      Int_t hitID = (fAPVKey << 8) | chNo ;
      Float_t data = 0;
      Float_t avgdata = 0;
      if(fMapping->GetPRadStripMapping(fAPVID, chNo)<16)
	{
	  data    = -(rawDataTS[chNo] - fPedestalOffsets[chNo] - commonModeOffsets_1stSet[timebin]) ;      
	  avgdata = -rawDataZS_1stSet[i1st] ;
	  i1st++;
	}
      else
	{
	  data    = -(rawDataTS[chNo] - fPedestalOffsets[chNo] - commonModeOffsets_2ndSet[timebin]) ;
	  avgdata = -rawDataZS_2ndSet[i2nd] ;
	  i2nd++;
	}

      // PERFORM ZERO SUPPRESSION
      
      if( fZeroSupCut > 0 ) {
	if( avgdata > (fZeroSupCut * fPedestalNoises[chNo]) ) {
	  if(!fListOfHitsClean[hitID]) {
	    GEMHit * hit = new GEMHit(hitID, fAPVID, chNo, fZeroSupCut,  fIsHitMaxOrTotalADCs) ;
	    fListOfHitsClean[hitID] = hit ;
	  }
	  fListOfHitsClean[hitID]->AddTimeBinADCs(timebin, data) ;
	}
      }
      // NO ZERO SUPPRESSION
      if( fSaveRawHit > 0) {
	if(!fListOfHits[hitID]) {
	  GEMHit * hit = new GEMHit(hitID, fAPVID, chNo, 0,  fIsHitMaxOrTotalADCs) ;
	  fListOfHits[hitID] = hit ;
	}
	fListOfHits[hitID]->AddTimeBinADCs(timebin, data) ;
      }
     
    }    
    firstdata = lastdata + 12 ;
    lastdata = firstdata + NCH ;
    rawDataTS.clear() ;
  }

  rawDataTS.clear() ;
  rawDataTS_1stSet.clear() ;
  rawDataTS_2ndSet.clear() ;
  dataTest_1stSet.clear() ;
  dataTest_2ndSet.clear() ;

  rawDataTS.clear() ;
  commonModeOffsets_1stSet.clear() ;
  commonModeOffsets_2ndSet.clear() ;
}

map < TString, list <GEMHit * > > GEMOnlineHitDecoder::GetListOfHitsFromPlanes() {
  map < Int_t, GEMHit * >::const_iterator hit_itr ;
  for(hit_itr = fListOfHits.begin(); hit_itr != fListOfHits.end(); ++hit_itr) { 
    GEMHit * hit = (* hit_itr).second ;
    TString planename = hit->GetPlane() ;
    fListOfHitsFromPlane[planename].push_back(hit) ;
    int adc = hit->GetHitADCs() ;
  }
  return fListOfHitsFromPlane ;
}

TH1F* GEMOnlineHitDecoder::GetHit(TString str)
{
  TH1F * h1;

  int NN = fMapping->GetNbOfPlane();

  int nbDetector = fMapping->GetNbOfDetectors();

  for(int i=0;i<nbDetector;i++)
  {
    TString detectorName = fMapping->GetDetectorFromID(i);
    list<TString> planeList = fMapping->GetPlaneListFromDetector(detectorName);
    int nbPlane = planeList.size();

    list<TString>::iterator it;
    for(it=planeList.begin();it!=planeList.end();++it)
    {
      if(*it == str)
      {
        TString hh = detectorName+"_"+(*it)+"_hit_distribution";
	int N = (int) (fMapping->GetPlaneSize(*it))/0.4;
	h1 = new TH1F(hh, hh, N, -fMapping->GetPlaneSize(*it)/2, fMapping->GetPlaneSize(*it)/2 );
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

map < TString, list <GEMHit * > > GEMOnlineHitDecoder::GetListOfHitsCleanFromPlanes() {
  map < Int_t, GEMHit * >::const_iterator hit_itr ;
  for(hit_itr = fListOfHitsClean.begin(); hit_itr != fListOfHitsClean.end(); ++hit_itr) { 
    GEMHit * hit = (* hit_itr).second ;
    TString planename = hit->GetPlane() ;
    fListOfHitsCleanFromPlane[planename].push_back(hit) ;
    Float_t adc = hit->GetHitADCs() ;
  }
  return fListOfHitsCleanFromPlane ;
}

TH1F* GEMOnlineHitDecoder::GetCleanHit(TString str)
{
  TH1F * h2;

  int NN = fMapping->GetNbOfPlane();

  int nbDetector = fMapping->GetNbOfDetectors();

  for(int i=0;i<nbDetector;i++)
  {
    TString detectorName = fMapping->GetDetectorFromID(i);
    list<TString> planeList = fMapping->GetPlaneListFromDetector(detectorName);
    int nbPlane = planeList.size();

    list<TString>::iterator it;
    for(it=planeList.begin();it!=planeList.end();++it)
    {
      if(*it == str)
      {
        TString hh = detectorName+"_clean_"+(*it)+"_hit_distribution";
	h2 = new TH1F(hh, hh, 2000, -fMapping->GetPlaneSize(*it)/2-100, fMapping->GetPlaneSize(*it)/2+100 );
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

Bool_t GEMOnlineHitDecoder::IsADCchannelActive() {
  Bool_t isADCchannelActive = kFALSE ;
  PRDMapping * mapping = PRDMapping::GetInstance() ;
  fActiveADCchannels = mapping->GetActiveADCchannels(fFECID) ;
  if ( find ( fActiveADCchannels.begin(), fActiveADCchannels.end(), fADCChannel ) != fActiveADCchannels.end() ) {
    //fAPVName           = mapping->GetAPVFromID(fAPVID);
    fAPVKey            = mapping->GetAPVNoFromID(fAPVID);
    fAPVHeaderLevel    = mapping->GetAPVHeaderLevelFromID(fAPVID);
    //fAPVIndexOnPlane   = mapping->GetAPVIndexOnPlane(fAPVID);
    //fAPVOrientation    = mapping->GetAPVOrientation(fAPVID);
    //fPlane             = mapping->GetPlaneFromAPVID(fAPVID); 
    //fPlaneID           = mapping->GetPlaneID(fPlane) ;
    //fDetector          = mapping->GetDetectorFromPlane(fPlane) ;
    //fDetectorID        = mapping->GetDetectorID(fDetector) ;
    //fDetectorType      = mapping->GetDetectorTypeFromDetector(fDetector) ;
    //fReadoutBoard      = mapping->GetReadoutBoardFromDetector(fDetector) ;
    //fPlaneSize         = mapping->GetPlaneSize(fPlane);
    //fNbOfAPVsOnPlane   = mapping->GetNbOfAPVsOnPlane(fPlane);
    isADCchannelActive = kTRUE ;
  }
  return isADCchannelActive ;
}

/*****************************************************************************
 * **                    Compute Clusters Information                     ** *
 *****************************************************************************/
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

void GEMOnlineHitDecoder::ComputeClusters()
{
  map < TString, list <GEMHit*> >::const_iterator  hitsFromPlane_itr ;

  for (hitsFromPlane_itr = fListOfHitsCleanFromPlane.begin(); hitsFromPlane_itr != fListOfHitsCleanFromPlane.end(); ++hitsFromPlane_itr) {
    TString plane =  (*hitsFromPlane_itr).first ;
    //cout<<"ComputeClusters:  "<<plane<<endl;
    list <GEMHit*> hitsFromPlane = (*hitsFromPlane_itr).second ; 
    hitsFromPlane.sort(CompareStripNo) ;
    Int_t listSize = hitsFromPlane.size() ;

    if (listSize < fMinClusterSize) {
      fIsGoodClusterEvent = kFALSE ;
      continue ;
    }

    Int_t previousStrip = -2 ;
    Int_t clusterNo = -1 ;
    map<Int_t, GEMCluster *> clustersMap ;
    list <GEMHit *>::const_iterator hit_itr ;

    for (hit_itr = hitsFromPlane.begin(); hit_itr != hitsFromPlane.end(); hit_itr++) {
      GEMHit * hit =  * hit_itr ; 
      Int_t currentStrip = hit->GetStripNo() ;

      // remove first 16 strips (apv index 0 on X side) and last 16 strips (apv index 10 on X side)
      if( plane.Contains("X") && ( (currentStrip<16) || (currentStrip > 1391) )  ) continue;

      Int_t apvIndexOnPlane = hit->GetAPVIndexOnPlane();
      if(currentStrip != (previousStrip + 1)) {
	clusterNo++ ;
      }
      if(!clustersMap[clusterNo]) {
	clustersMap[clusterNo] = new GEMCluster(fMinClusterSize, fMaxClusterSize, fIsClusterMaxOrTotalADCs) ;
	clustersMap[clusterNo]->SetNbAPVsFromPlane(hit->GetNbAPVsFromPlane());
	clustersMap[clusterNo]->SetAPVIndexOnPlane(hit->GetAPVIndexOnPlane());
	clustersMap[clusterNo]->SetPlaneSize(hit->GetPlaneSize());
	clustersMap[clusterNo]->SetPlane(hit->GetPlane());
      }
      clustersMap[clusterNo]->AddHit(hit) ;
      previousStrip = currentStrip;
    }

    map<Int_t, GEMCluster *>::const_iterator  cluster_itr ;
    for (cluster_itr = clustersMap.begin(); cluster_itr != clustersMap.end(); cluster_itr++) {
      GEMCluster * cluster = ( * cluster_itr ).second ;
      if (!cluster->IsGoodCluster()) {
	delete cluster ;
	continue ;
      }
      cluster->ComputeClusterPosition() ;
      fListOfClustersCleanFromPlane[plane].push_back(cluster) ;
    }

    fListOfClustersCleanFromPlane[plane].sort(CompareClusterADCs) ;
    hitsFromPlane.clear() ;
    clustersMap.clear() ;
  }

}

TH1F* GEMOnlineHitDecoder::GetCluster(TString str)
{
  TH1F * hc1;

  int NN = fMapping->GetNbOfPlane();

  int nbDetector = fMapping->GetNbOfDetectors();

  for(int i=0;i<nbDetector;i++)
  {
    TString detectorName = fMapping->GetDetectorFromID(i);
    list<TString> planeList = fMapping->GetPlaneListFromDetector(detectorName);
    int nbPlane = planeList.size();

    list<TString>::iterator it;
    for(it=planeList.begin();it!=planeList.end();++it)
    {
      if(*it == str)
      {
        TString hh = detectorName+"_"+(*it)+"_Cluster_Distribution";
	hc1 = new TH1F(hh, hh, 2000, -fMapping->GetPlaneSize(*it)/2-100, fMapping->GetPlaneSize(*it)/2+100 );
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

void GEMOnlineHitDecoder::DeleteClustersInPlaneMap() {
  map < TString, list <GEMCluster *> >::const_iterator itr ;
  for (itr = fListOfClustersCleanFromPlane.begin(); itr != fListOfClustersCleanFromPlane.end(); itr++) {
    list <GEMCluster *> listOfCluster = (*itr).second ;
    listOfCluster.clear() ;
  }
  fListOfClustersCleanFromPlane.clear() ;
}

void GEMOnlineHitDecoder::GetClusterHyCalCutMode(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersCleanFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersCleanFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersCleanFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersCleanFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  int s2 = cluster_y1.size();  
  int nbCluster1 = (s1<s2)?s1:s2;
  s1 = cluster_x2.size();
  s2 = cluster_y2.size();
  int nbCluster2 = (s1<s2)?s1:s2;

  /*
   * Do the coordinate convert here
   * Convert GEM coordinate to HyCal coordinate
   *
   * **************************************
   *    Move GEM coordinate to Beam Hole center
   *    overlapping area: hole diameter (44mm)
   *    origin shift: 550.4/2 - (44-pitch)/2 = 253.2
   *
   *    gem1: x = x-253.2; y = y
   *    gem2: x = -x+253.2; y=-y
   *    
   *    right-hand coordinate, HyCal Y axis
   *    must pointing downward
   *
   *    beam downstream is z axis direction
   * **************************************
   *
   */
  double O_Transfer = 253.2;
  double OverlapLength = 44;
  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm


   // cutting edge
  float edge1 = 0;
  float edge2 = 0;

  // offset from data
  double xoffset = -0.3618;
  double yoffset = 0.1792;

  //the above offsets are from the projection on GEM2
  //real GEM1 offsets should be projected back
  xoffset = xoffset*z_gem1/z_gem2;
  yoffset = yoffset*z_gem1/z_gem2;

  if(nbCluster1>0)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    for(int i = 0;i<nbCluster1;i++)
    {
      if(((*itx)->GetClusterPosition() -xoffset -O_Transfer) <= edge1) 
      // remove overlapping area on GEM1, use the corresponding area on GEM2
      // do not use edge as the cut line. choose some arbitrary line and find the corresponding line on GEM2
      {
        float c_x = (*itx)->GetClusterADCs();
	float c_y = (*ity)->GetClusterADCs();
	float x = (*itx++)->GetClusterPosition() -O_Transfer -xoffset;
	float y = (*ity++)->GetClusterPosition() -yoffset; 
	gem1.push_back(GEMClusterStruct(x, y, c_x, c_y));
      }
    }
  }

   if(nbCluster2>0)
  {
    list<GEMCluster*>::iterator itx2 = cluster_x2.begin();
    list<GEMCluster*>::iterator ity2 = cluster_y2.begin();
    for(int i = 0;i<nbCluster2;i++)
    {
      if( ( O_Transfer - ((*itx2)->GetClusterPosition())  ) > edge2)
      {
        float c_x = (*itx2)->GetClusterADCs();
	float c_y = (*ity2)->GetClusterADCs();
	float x = O_Transfer - ((*itx2++)->GetClusterPosition()); 
	float y =  -(*ity2++)->GetClusterPosition(); 
	gem2.push_back(GEMClusterStruct(x, y, c_x, c_y));
      }
    }
  }

}

void GEMOnlineHitDecoder::GetClusterHyCalPlusMode(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersCleanFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersCleanFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersCleanFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersCleanFromPlane["pRadGEM2Y"];

  /*
   * Do the coordinate convert here
   * Convert GEM coordinate to HyCal coordinate
   *
   * **************************************
   *    Move GEM coordinate to HyCal Hole center
   *    overlapping area: hole diameter (44mm)
   *    origin shift: 550.4/2 - (44-pitch)/2 = 253.2
   *
   *    gem1: x = x-253.2; y = y
   *    gem2: x = -x+253.2; y=-y
   *    
   *    right-hand coordinate, HyCal Y axis
   *    must pointing downward
   *
   *    beam downstream is z axis direction
   * **************************************
   *
   */
  double O_Transfer = 253.2;
  double OverlapLength = 44;
  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm

  // offset from data
  double xoffset = -0.3618;
  double yoffset = 0.1792;

  //the above offsets are from the projection on GEM2
  //real GEM1 offsets should be projected back
  xoffset = xoffset*z_gem1/z_gem2;
  yoffset = yoffset*z_gem1/z_gem2;

  for( auto &i : cluster_x1)
  {
      for(auto &j : cluster_y1)
      {
          float c_x = i->GetClusterADCs();
	  float c_y = j->GetClusterADCs();
	  float x = i->GetClusterPosition() - O_Transfer-xoffset;
	  float y = j->GetClusterPosition() - yoffset;
	  gem1.push_back(GEMClusterStruct(x,y,c_x,c_y));
      }
  }
  for( auto &i : cluster_x2)
  {
      for(auto &j : cluster_y2)
      {
          float c_x = i->GetClusterADCs();
	  float c_y = j->GetClusterADCs();
	  float x = O_Transfer - i->GetClusterPosition();
	  float y = - j->GetClusterPosition() ;
	  gem2.push_back(GEMClusterStruct(x,y,c_x,c_y));
      }
  }

}

void GEMOnlineHitDecoder::GetClusterBeamLine(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersCleanFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersCleanFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersCleanFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersCleanFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  int s2 = cluster_y1.size();  
  int nbCluster1 = (s1<s2)?s1:s2;
  s1 = cluster_x2.size();
  s2 = cluster_y2.size();
  int nbCluster2 = (s1<s2)?s1:s2;

  /*
   * Do the coordinate convert here
   * Convert GEM coordinate to HyCal coordinate
   *
   * **************************************
   *    Move GEM coordinate to Beam Hole center
   *    overlapping area: hole diameter (44mm)
   *    origin shift: 550.4/2 - (44-pitch)/2 = 253.2
   *
   *    gem1: x = x-253.2; y = y
   *    gem2: x = -x+253.2; y=-y
   *    
   *    right-hand coordinate, HyCal Y axis
   *    must pointing downward
   *
   *    beam downstream is z axis direction
   * **************************************
   *
   */
  double O_Transfer = 253.2;
  double OverlapLength = 44;
 
  // cutting edge
  float edge1 = 0;
  float edge2 = 0;
  double z_gem1 = 5300; //mm
  double z_gem2 = 5260; //mm

  double xoffset_gem = -0.3618;
  double yoffset_gem = 0.1792;

  //the above offsets are from the projection on GEM2
  //real GEM1 offsets should be projected back
  xoffset_gem = xoffset_gem*z_gem1/z_gem2;
  yoffset_gem = yoffset_gem*z_gem1/z_gem2;


  double xoffset_beam = 1.631;
  double yoffset_beam = 0.366;

  if(nbCluster1>0)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    for(int i = 0;i<nbCluster1;i++)
    {
      /*
       * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       * overlapping area: according to frame design, overlapping area is in fact 44mm
       * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       */
      if(((*itx)->GetClusterPosition() -xoffset_gem -O_Transfer) <= edge1) 
      // remove overlapping area on GEM1, use the corresponding area on GEM2
      // do not use edge as the cut line. choose some arbitrary line and find the corresponding line on GEM2
      {
        float c_x = (*itx)->GetClusterADCs();
	float c_y = (*ity)->GetClusterADCs();
	float x = ((*itx++)->GetClusterPosition()) -O_Transfer - xoffset_gem - xoffset_beam;
	float y = (*ity++)->GetClusterPosition() -yoffset_gem - yoffset_beam;
	gem1.push_back(GEMClusterStruct(x, y, c_x, c_y));
      }
    }
  }

   if(nbCluster2>0)
  {
    list<GEMCluster*>::iterator itx2 = cluster_x2.begin();
    list<GEMCluster*>::iterator ity2 = cluster_y2.begin();
    for(int i = 0;i<nbCluster2;i++)
    {
      if( ( O_Transfer - ((*itx2)->GetClusterPosition())  ) < edge2) continue;
      {
        float c_x = (*itx2)->GetClusterADCs();
	float c_y = (*ity2)->GetClusterADCs();
	float x = O_Transfer - ((*itx2++)->GetClusterPosition()) - xoffset_beam;
	float y =  -(*ity2++)->GetClusterPosition()  - yoffset_beam;
	gem2.push_back(GEMClusterStruct(x, y, c_x, c_y));
      }
    }
  }

}

void GEMOnlineHitDecoder::GetClusterGEM(vector<GEMClusterStruct> &gem1, vector<GEMClusterStruct> &gem2)
{
  list<GEMCluster*> cluster_x1 = fListOfClustersCleanFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersCleanFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersCleanFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersCleanFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  int s2 = cluster_y1.size();
  int nbCluster1 = (s1<s2)?s1:s2;
  s1 = cluster_x2.size();
  s2 = cluster_y2.size();
  int nbCluster2 = (s1<s2)?s1:s2;

  if(nbCluster1>0)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    for(int i = 0;i<nbCluster1;i++)
    {
      float c_x = (*itx)->GetClusterADCs();
      float c_y = (*ity)->GetClusterADCs();
      float x = (*itx++)->GetClusterPosition();
      float y = (*ity++)->GetClusterPosition(); 
      gem1.push_back(GEMClusterStruct(x, y, c_x, c_y));
    }
  }

   if(nbCluster2>0)
  {
    list<GEMCluster*>::iterator itx2 = cluster_x2.begin();
    list<GEMCluster*>::iterator ity2 = cluster_y2.begin();
    for(int i = 0;i<nbCluster2;i++)
    {
      float c_x = (*itx2)->GetClusterADCs();
      float c_y = (*ity2)->GetClusterADCs();
      float x =  (*itx2++)->GetClusterPosition() ; 
      float y =  (*ity2++)->GetClusterPosition() ; 
      gem2.push_back(GEMClusterStruct(x, y, c_x, c_y));
    }
  }

}

void GEMOnlineHitDecoder::FillHistos(TH1F* hNbClusterPerPlaneX[], TH1F* hNbClusterPerPlaneY[], TH1F* hClusterDistX[], TH1F* hClusterDistY[])
{
  list<GEMCluster*> cluster_x1 = fListOfClustersCleanFromPlane["pRadGEM1X"];
  list<GEMCluster*> cluster_y1 = fListOfClustersCleanFromPlane["pRadGEM1Y"];
  list<GEMCluster*> cluster_x2 = fListOfClustersCleanFromPlane["pRadGEM2X"];
  list<GEMCluster*> cluster_y2 = fListOfClustersCleanFromPlane["pRadGEM2Y"];

  int s1 = cluster_x1.size();
  hNbClusterPerPlaneX[0]->Fill(s1);
  int s2 = cluster_y1.size();
  hNbClusterPerPlaneY[0]->Fill(s2);
  for(int i=0;i<s1;i++)
  {
    list<GEMCluster*>::iterator itx = cluster_x1.begin();
    hClusterDistX[0] -> Fill( ( *itx++)->GetClusterPosition() );
  }
  for(int i=0;i<s2;i++)
  {
    list<GEMCluster*>::iterator ity = cluster_y1.begin();
    hClusterDistY[0]->Fill( (*ity++)->GetClusterPosition() );
  }

  s1 = cluster_x2.size();
  hNbClusterPerPlaneX[1] -> Fill(s1);
  s2 = cluster_y2.size();
  hNbClusterPerPlaneY[1] -> Fill(s2);
  for(int i=0;i<s1;i++)
  {
    list<GEMCluster*>::iterator itx = cluster_x2.begin();
    hClusterDistX[1] -> Fill( ( *itx++)->GetClusterPosition() );
  }
  for(int i=0;i<s2;i++)
  {
    list<GEMCluster*>::iterator ity = cluster_y2.begin();
    hClusterDistY[1]->Fill( (*ity++)->GetClusterPosition() );
  }
}
