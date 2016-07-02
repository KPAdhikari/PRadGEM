#include "GEMHitDecoder.h" 
//ClassImp (GEMHitDecoder);

static string trim(const string &str, const string &w = " \t\n\r")
{

    const auto strBegin = str.find_first_not_of(w);
    if (strBegin == string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(w);

    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


//==================================================================================
GEMHitDecoder::GEMHitDecoder(string str) {
 //   printf("Enter GEMHitDecoder \n");
  file_name = str;
  NCH = 128;
  fNbOfChannels = NCH; 
  //fNWords = 100 ;
  fEventNb = -1 ;
  //fPacketSize = 4000 ;
  //fIsFirstEvent = kTRUE;
  //fIsNewPacket = kFALSE ;
  //
  //from gem configure file
  //fZeroSupCut = zeroSup ;
  fZeroSupCut = 5;
  //fZeroSupCut = 0;
  fRunType = "HITS" ;
  fCurrentFEC = -1 ;
  fIsNewFECdataFlag = kFALSE ;
  fStartData = 0 ;
  fIsHitMaxOrTotalADCs = "signalPeak" ;

  fAPVID = 16 ;
  fAPVName           = "defaultAPV" ;
  fAPVKey            = 0 ;
  fAPVHeaderLevel    = 1200 ;
  fAPVIndexOnPlane   = 0;
  fAPVOrientation    = 0;
  fPlane             = "PRADGEMX" ;
  fDetector          = "PRADGEM" ;
  fDetectorID        = 0 ;
  fDetectorType      = "PRADGEM" ;
  fReadoutBoard      = "CARTESIAN" ;
  fPlaneID           = 0 ;
  fPlaneSize         = 0.0 ;
  fNbOfAPVsOnPlane   = 1 ;

  //xb
  fMapping = PRDMapping::GetInstance();
  FECs.clear();
  mSrsSingleEvent.clear();
  FECs = fMapping->GetBankIDSet();
  cout<<"FECs Found: ";
  set<int>::iterator it;
  for(it=FECs.begin(); it!=FECs.end();++it)
  {
    cout<<(*it)<<"  ";
  }
  cout<<endl;
  // to be read from configure file
  pedestal_file = "./pedestal.root";
  //TFile *f = new TFile(pedestal_file.c_str() );
  ped = new GEMPedestal(pedestal_file);
  ped->LoadPedestal();

  //  printf("Exit  GEMHitDecoder \n");
}
/*
//==================================================================================
GEMHitDecoder::GEMHitDecoder() {
  //  printf("Enter GEMHitDecoder \n");
  NCH = 128;
  fNbOfChannels = NCH; 
  //fNWords = 100 ;
  fEventNb = -1 ;
  //fPacketSize = 4000 ;
  //fIsFirstEvent = kTRUE;
  //fIsNewPacket = kFALSE ;
  fZeroSupCut = 0 ;
  fRunType = "HITS" ;
  fCurrentFEC = -1 ;
  fIsNewFECdataFlag = kFALSE ;
  fStartData = 0 ;
  fIsHitMaxOrTotalADCs = "signalPeak" ;

  fAPVID = 16 ;
  fAPVName           = "defaultAPV" ;
  fAPVKey            = 0 ;
  fAPVHeaderLevel    = 1200 ;
  fAPVIndexOnPlane   = 0;
  fAPVOrientation    = 0;
  fPlane             = "PRADGEMX" ;
  fDetector          = "PRADGEM" ;
  fDetectorID        = 0 ;
  fDetectorType      = "PRADGEM" ;
  fReadoutBoard      = "CARTESIAN" ;
  fPlaneID           = 0 ;
  fPlaneSize         = 0.0 ;
  fNbOfAPVsOnPlane   = 1 ;
  //  printf("Exit  GEMHitDecoder %d \n", fZeroSupCut);
}
*/
//==================================================================================
GEMHitDecoder::~GEMHitDecoder() {
  Clear() ;
  fListOfHits.clear() ;
  fListOfHitsFromPlane.clear() ;
  fActiveADCchannels.clear() ;
}

//==================================================================================
void GEMHitDecoder::Clear() {
  fRawData16bits.clear();
  fPedestalNoises.clear(); 
  fPedestalOffsets.clear();
}

//===================================================================
Bool_t GEMHitDecoder::IsADCchannelActive() {
  Bool_t isADCchannelActive = kFALSE ;
  PRDMapping * mapping = PRDMapping::GetInstance() ;
  fActiveADCchannels = mapping->GetActiveADCchannels(fFECID) ;
  if ( find ( fActiveADCchannels.begin(), fActiveADCchannels.end(), fADCChannel ) != fActiveADCchannels.end() ) {
    fAPVName           = mapping->GetAPVFromID(fAPVID);
    fAPVKey            = mapping->GetAPVNoFromID(fAPVID);
    fAPVHeaderLevel    = mapping->GetAPVHeaderLevelFromID(fAPVID);
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
    isADCchannelActive = kTRUE ;
  }
  return isADCchannelActive ;
}

//=========================================================================================
Int_t GEMHitDecoder::ProcessAllEvents(Int_t evtID)
{
  cout<<"Enter GEMHitDecoder::ProcessAllEvents():  "<<endl;
  int entry = 0;
  vector<int> event_vector;

  try
  {
    evioFileChannel chan(file_name.c_str(), "r");
    chan.open();
    while(chan.read())
    {
      //xb
      Clear();
      fListOfHits.clear();
      fListOfHitsFromPlane.clear();

      //vector<int> event_vector;
      event_vector.clear();
      event_vector.resize(0);

      evioDOMTree event(chan);
      evioDOMNodeListP eventList = event.getNodeList( isLeaf() );
      cout<<"total number of all banks: "<<eventList->size()<<endl;
      
      evioDOMNodeList::iterator iter;
      for(iter=eventList->begin(); iter!=eventList->end(); ++iter)
      {
        if( FECs.find( ((*iter)->tag)-9 ) != FECs.end() )
	{
	  vector<uint32_t> *vec = (*iter)->getVector<uint32_t>();
	  if(vec!=NULL && vec->size()>0)
	  {
	    event_vector.reserve( event_vector.size() + vec->size() );
	    event_vector.insert(event_vector.end(), vec->begin(), vec->end() );
	  }
	  else
	  {
	    cout<<"Found NULL contents in fec ..."<<endl;
	  }
	}
      }

      ++entry;
      cout<<"entry: "<<entry<<endl;
      cout<<"SRS Event Size[uint_32]:"<<event_vector.size()<<endl;

      //cout<<"evtID:  "<<evtID<<endl;
      if(entry >= evtID && evtID!=-1) break;
      
      if(event_vector.size() == 0) continue;

      GEMRawDecoder raw_decoder(event_vector);

      //cout<<"reached here..."<<endl;
      EventHandler( raw_decoder.GetDecoded(), ped);
      //cout<<"reached EventHandler here..."<<endl;

      ShowHits();
    }
  } catch (evioException e) {
    cerr<<endl<<e.toString()<<endl<<endl;
    exit(EXIT_FAILURE);
  }
  return entry;
}

//===================================================================
void GEMHitDecoder::EventHandler( map<int, map<int, vector<int> > > srsSingleEvent , GEMPedestal * ped )
{
  //cout<<"GEMHitDecoder::EventHandler begin..."<<endl;

  mSrsSingleEvent = srsSingleEvent;
  //cout<<"single event size:  "<<mSrsSingleEvent.size()<<endl;
  map<int, map<int, vector<int> > >::iterator it;
  for(it=mSrsSingleEvent.begin(); it!=mSrsSingleEvent.end(); ++it)
  {
    //cout<<"xb: APVs found in current FEC: "<<mSrsSingleEvent.size()<<endl;
    fFECID = it->first;
    map<int, vector<int> > fec_data = it->second;
    map<int, vector<int> >::iterator itt;
    //cout<<"xb: APVs found in current FEC: "<<fec_data.size()<<endl;
    for(itt=fec_data.begin(); itt!=fec_data.end(); ++itt)
    {
      fADCChannel = itt->first;
      fAPVID = fFECID << 4 | fADCChannel;
      if(IsADCchannelActive())
      {
        fRawData16bits.clear();
	fRawData16bits = itt->second;

	//Get Apv pedestals ...
	fPedestalNoises.resize(NCH);
	fPedestalOffsets.resize(NCH);
	fPedestalNoises = ped->GetAPVNoises(fAPVID);
	fPedestalOffsets = ped->GetAPVOffsets(fAPVID);
	fPedestalNoises_1stSet.clear();
	fPedestalOffsets_1stSet.clear();
	fPedestalNoises_2ndSet.clear();
	fPedestalOffsets_2ndSet.clear();

	assert( fRawData16bits.size() > 0 );
	
	// Decode APV25 DATA INTO HITS
	fAPVStatus = fMapping->GetAPVstatus(fAPVID);
	string apv_string = fAPVStatus.Data();
	apv_string = trim(apv_string);
	//cout<<fAPVStatus<<endl;
	if( apv_string != "normal" )
	{
	  //cout<<"xb: !normal"<<endl;
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
	  //cout<<"xb: normal"<<endl;
          APVEventDecoder();
	}
	//fPedestalNoises.clear();
	//fPedestalOffsets.clear();
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

}

//=====================================================
void GEMHitDecoder::APVEventDecoder() {
  //printf(" Enter  GEMHitDecoder::APVEventDecoder()\n") ;

  Int_t idata = 0, firstdata = 0, lastdata = 0, nbOfTimeBin = 3;
  Int_t size = fRawData16bits.size() ;
  //cout<<"   APVEventDecoder:: apv size: "<<size<<endl;
  Float_t commonMode = 0 ;

  vector<Float_t> rawDataTS, rawDataZS;
  vector<Float_t> threshold ;
  vector<Float_t> dataTest, commonModeOffsets ;
  rawDataTS.clear();

  // COMPUTE APV25 COMMON MODE CORRECTION
  //if (fIsNewFECdataFlag) {
  if(1){
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
    //fIsNewFECdataFlag = kFALSE ;
  }
  //xb check raw data
  /*
  TH1F *hh = new TH1F("hh", "hh", size, 0, size-1);
  hh->SetStats(0);
  for(int j=0;j<size;j++)
  {
    hh->SetBinContent(j+1, fRawData16bits[j] ); 
  }
  */

  //===========================================================================================================
  // BLOCK WHERE COMMON MODE CORRECTION FOR ALL APV25 TIME BINS IS COMPUTED
  //===========================================================================================================
  firstdata = fStartData ;
  lastdata  = firstdata  + NCH ;
  //cout<<"total size:  "<<size<<endl;
  //cout<<"first data: "<<firstdata<< "  last data: "<<lastdata<<endl;

  for(Int_t timebin = 0; timebin < nbOfTimeBin; timebin++) {

    // EXTRACT APV25 DATA FOR A GIVEN TIME BIN
    if( lastdata > size ) 
    {
      cout<<"###Warning:  GEMHitDecoder: APV Raw Data out of Range"<<endl;
      cout<<"           please adjust your header level in mapping file"<<endl;
      cout<<"           increase recommended..."<<endl;
    }
    rawDataTS.insert(rawDataTS.end(), &fRawData16bits[firstdata], &fRawData16bits[lastdata]);
    assert(rawDataTS.size() == 128 );

    //xb check raw data
    /*
    TCanvas *c = new TCanvas("c", "c", 10, 10, 1200, 400);
    c->Divide(2,1);
    TH1F *h = new TH1F("h", "h", 128, 0.5, 127.5);
    h->GetXaxis()->SetRange(0,3500);
    h->SetStats(0);
    for(int i=0;i<NCH;i++)
    {
      h->SetBinContent(i+1, fRawData16bits[firstdata+i]);
    }
    */

    // PERFORM APV25 PEDESTAL OFFSET CORRECTION  FOR A GIVEN TIME BIN
    std::transform (rawDataTS.begin(), rawDataTS.end(), fPedestalOffsets.begin(), rawDataTS.begin(), std::minus<Float_t>());

    // COMPUTE THRESOLD (= 20 PEDESTAL NOISE) FOR COMMON MODE DATA SELECTION 
    threshold.resize(rawDataTS.size()) ;
    std::transform (fPedestalNoises.begin(), fPedestalNoises.end(), threshold.begin(), std::bind1st(std::multiplies<Float_t>(), 20.0));

    dataTest.clear();
    dataTest.insert(dataTest.begin(), rawDataTS.begin(), rawDataTS.end() );
    int nbHits = 0;
    assert( rawDataTS.size() == 128 );
    for(int i=0; i<128; i++)
    {
      if( (dataTest[i]-threshold[i])>=0 )
      {
        dataTest[i] = 0;
	nbHits++;
      }
    }
    assert( nbHits <= 128);

    // COMPUTE COMMON MODE FOR A GIVEN APV AND TIME BIN
    commonMode = std::accumulate(dataTest.begin(), dataTest.end(), 0.0) / (NCH - nbHits) ;
    commonModeOffsets.push_back(commonMode) ;
    //xb check here
    
    //cout<<"signal strips:  "<<nbHits<<endl;
    //cout<<"TS:  "<< timebin<<"  commonMode:  "<<commonMode<<endl;
    /*
    c->cd(1);
    hh->Draw();
    c->cd(2);
    h->Draw();
    */

    // PERFORM COMMON MODE CORRECTION FOR A GIVEN TIME BIN
    std::transform (rawDataTS.begin(), rawDataTS.end(), rawDataTS.begin(), std::bind2nd(std::minus<Float_t>(), commonMode));

    //  ADC SUM OVER ALL TIME BINS USE AS THE TEST CRITERIA FOR ZERO SUPPRESSION
    if(timebin == 0)  rawDataZS.resize(rawDataTS.size()) ;
    std::transform (rawDataZS.begin(), rawDataZS.end(), rawDataTS.begin(), rawDataZS.begin(), std::plus<Float_t>());

    // PROCEED TO NEXT TIME BIN
    firstdata = lastdata + 12 ;
    lastdata = firstdata + NCH ;

    // CLEAR EVERYTHING
    rawDataTS.clear() ;
    threshold.clear() ;
    dataTest.clear() ;
    /*
    c->Update();
    getchar();
    c->Delete();
    h->Delete();
    */
  }
  //getchar();


  // ADC AVERAGE OVER ALL TIME BINS USE AS THE TEST CRITERIA FOR ZERO SUPPRESSION
  std::transform (rawDataZS.begin(), rawDataZS.end(), rawDataZS.begin(), std::bind2nd(std::divides<Float_t>(), nbOfTimeBin)) ;

  //=============================================================================================================================
  // BLOCK WHERE ZERO SUPPRESSION IS COMPUTED
  //=============================================================================================================================
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
	  if(!fListOfHits[hitID]) 
	  {
	    GEMHit * hit = new GEMHit(hitID, fAPVID, chNo, fZeroSupCut,  fIsHitMaxOrTotalADCs) ;
	    fListOfHits[hitID] = hit ;
	  }
	  fListOfHits[hitID]->AddTimeBinADCs(timebin, data) ;
	}
      }

      // NO ZERO SUPPRESSION
      else 
      {
	Float_t adc = data ;
	if( avgdata < (5 * fPedestalNoises[chNo]) ) adc = avgdata ;
	if(!fListOfHits[hitID]) 
	{
	  GEMHit * hit = new GEMHit(hitID, fAPVID, chNo, fZeroSupCut,  fIsHitMaxOrTotalADCs) ;
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

//===================================================================
map < TString, list <GEMHit * > > GEMHitDecoder::GetListOfHitsFromPlanes() {
  // printf(" Enter  GEMHitDecoder::GetListOfHitsFromPlanes()\n") ;
  map < Int_t, GEMHit * >::const_iterator hit_itr ;
  for(hit_itr = fListOfHits.begin(); hit_itr != fListOfHits.end(); ++hit_itr) { 
    GEMHit * hit = (* hit_itr).second ;
    TString planename = hit->GetPlane() ;
    fListOfHitsFromPlane[planename].push_back(hit) ;
    Float_t adc = hit->GetHitADCs() ;
  }
  return fListOfHitsFromPlane ;
}

//======================================================================
void GEMHitDecoder::ShowHits()
{
  //ProcessAllEvents();
  GetListOfHitsFromPlanes();
  Int_t NN = fMapping->GetNbOfPlane();
  //cout<<"total number of planes:  "<<NN<<endl;
  TH1F* h[NN];

  Int_t nbDetector = fMapping->GetNbOfDetectors();
  //cout<<"nb of  detectors:  "<<nbDetector<<endl;
  for(int i=0;i<nbDetector;i++)
  {
    TString detectorName = fMapping->GetDetectorFromID(i);
    list<TString> planeList = fMapping->GetPlaneListFromDetector(detectorName);
    Int_t nbPlane = planeList.size();
    //cout<< "nb of plane: "<<nbPlane<<endl;
    list<TString>::iterator it;
    int kk = 0;
    for( it=planeList.begin();it!=planeList.end();++it)
    {
      //cout<<*it<<endl;
      //stringstream out;
      //out<<
      TString hh = detectorName+"_"+(*it);
      h[ i*nbPlane + kk ] = new TH1F(hh, hh, 2000, -fMapping->GetPlaneSize(*it)/2-100, fMapping->GetPlaneSize(*it)/2+100 );
      list< GEMHit* > hitList = fListOfHitsFromPlane[ *it  ];
      //cout<<"nb of hits on plane: "<<hitList.size()<<endl;
      list< GEMHit* >::iterator hit_it;
      for(hit_it=hitList.begin(); hit_it!=hitList.end();++hit_it)
      {
        Float_t pos = (*hit_it) -> GetStripPosition();
	Float_t adc = (*hit_it) -> GetHitADCs(); 
	h[i*nbPlane + kk] -> Fill(pos, adc);
      }
      kk++;
    }
  }
  TCanvas *cc = new TCanvas("cc", "cc", 10, 10, 1000, 1000);
  cc->Divide(1,4);
  for(int ii = 0; ii<NN; ii++)
  {
    cc->cd(ii+1);
    h[ii]->SetStats(0);
    //h[ii]->GetYaxis()->SetRange(0, 1000);
    h[ii]->SetMinimum(-500);
    h[ii]->SetMaximum(1200);
    h[ii]->Draw();
  }
  cc->Update();
  cc->Print("hit.png");
  getchar();
  for(int jj=0;jj<NN;jj++)
  {
    h[jj]->Delete();
  }
  //c->Delete();

  cout<<"number of hits: "<<fListOfHits.size()<<endl;
}


//=====================================================
void GEMHitDecoder::APVEventSplitChannelsDecoder() {
  //printf(" Enter  GEMHitDecoder::APVEventDecoderSpecialChannel()\n") ;

  Int_t idata = 0, firstdata = 0, lastdata = 0, nbOfTimeBin = 3;
  Int_t size = fRawData16bits.size() ;
  Int_t size1 = fPedestalNoises_1stSet.size();
  Int_t size2 = fPedestalNoises_2ndSet.size();

  //cout<<"   APVEventDecoder:: apv size: "<<size<<endl;
  Float_t commonMode = 0 ;

  vector<Float_t> rawDataTS, rawDataTS_1stSet, rawDataTS_2ndSet; 
  vector<Float_t> rawDataZS_1stSet, rawDataZS_2ndSet;
  vector<Float_t> threshold_1stSet, threshold_2ndSet ;
  vector<Float_t> dataTest_1stSet, dataTest_2ndSet;
  vector<Float_t> commonModeOffsets_1stSet, commonModeOffsets_2ndSet ;

  rawDataTS.clear(); 
  rawDataTS_1stSet.clear();
  rawDataTS_2ndSet.clear();

  //rawDataZS_1stSet.resize(rawDataTS_1stSet.size()) ;
  //rawDataZS_2ndSet.resize(rawDataTS_2ndSet.size()) ;
  threshold_1stSet.clear();
  threshold_2ndSet.clear();
  dataTest_1stSet.clear();
  dataTest_2ndSet.clear();
  commonModeOffsets_1stSet.clear();
  commonModeOffsets_2ndSet.clear();

  // COMPUTE APV25 COMMON MODE CORRECTION
  //if (fIsNewFECdataFlag) {
  //if(1){
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
    //fIsNewFECdataFlag = kFALSE ;
  //}
  //xb check raw data
  /*
  TH1F *hh = new TH1F("hh", "hh", size, 0, size-1);
  hh->SetStats(0);
  for(int j=0;j<size;j++)
  {
    hh->SetBinContent(j+1, fRawData16bits[j] ); 
  }
  */

  //===========================================================================================================
  // BLOCK WHERE COMMON MODE CORRECTION FOR ALL APV25 TIME BINS IS COMPUTED
  //===========================================================================================================
  firstdata = fStartData ;
  lastdata  = firstdata  + NCH ;
  //cout<<"total size:  "<<size<<endl;
  //cout<<"first data: "<<firstdata<< "  last data: "<<lastdata<<endl;

  for(Int_t timebin = 0; timebin < nbOfTimeBin; timebin++) {

    threshold_1stSet.resize(size1) ;
    threshold_2ndSet.resize(size2) ; 
    dataTest_1stSet.resize(size1) ;
    dataTest_2ndSet.resize(size2) ;
    rawDataZS_1stSet.resize(size1) ;
    rawDataZS_2ndSet.resize(size2) ;

    //=======================================================
    // EXTRACT APV25 DATA FOR A GIVEN TIME BIN
    //=======================================================
    if( lastdata > size ) 
    {
      cout<<"###ERROR:  GEMHitDecoder: APV Raw Data out of Range"<<endl;
      cout<<"           please adjust your header level in mapping file"<<endl;
      cout<<"           increase recommended..."<<endl;
    }
    assert( lastdata < size);
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

    //xb check raw data
    /*
    TCanvas *c = new TCanvas("c", "c", 10, 10, 1200, 400);
    c->Divide(2,1);
    TH1F *h = new TH1F("h", "h", 128, 0.5, 127.5);
    h->GetXaxis()->SetRange(0,3500);
    h->SetStats(0);
    for(int i=0;i<NCH;i++)
    {
      h->SetBinContent(i+1, fRawData16bits[firstdata+i]);
    }
    */

    //===============================================================
    //                     1ST SET OF CHANNELS
    //===============================================================

    // PERFORM APV25 PEDESTAL OFFSET CORRECTION  FOR A GIVEN TIME BIN
    std::transform (rawDataTS_1stSet.begin(), rawDataTS_1stSet.end(), fPedestalOffsets_1stSet.begin(), rawDataTS_1stSet.begin(), std::minus<Float_t>());

    // COMPUTE THRESOLD (= 20 PEDESTAL NOISE) FOR COMMON MODE DATA SELECTION 
    //threshold.resize(rawDataTS.size()) ;
    std::transform (fPedestalNoises_1stSet.begin(), fPedestalNoises_1stSet.end(), threshold_1stSet.begin(), std::bind1st(std::multiplies<Float_t>(), 50.0));

    //if(timebin == 0) dataTest.resize(rawDataTS.size()) ;
    dataTest_1stSet.clear();
    dataTest_1stSet.insert(dataTest_1stSet.begin(), rawDataTS_1stSet.begin(), rawDataTS_1stSet.end());
    assert( rawDataTS_1stSet.size() == 16 );
    int nbHits_1stSet = 0;
    for(int i=0;i<16;i++)
    {
      if( (dataTest_1stSet[i] - threshold_1stSet[i]) >= 0)
      {
        dataTest_1stSet[i] = 0;
	nbHits_1stSet++;
      }
    }
    assert(nbHits_1stSet <= 16 );

    // COMPUTE COMMON MODE FOR A GIVEN APV AND TIME BIN
    commonMode = std::accumulate(dataTest_1stSet.begin(),dataTest_1stSet.end(), 0.0) / (16 - nbHits_1stSet) ;
    commonModeOffsets_1stSet.push_back(commonMode) ;
    //xb check here
    
    //cout<<"signal strips:  "<<nbHits<<endl;
    //cout<<"TS:  "<< timebin<<"  commonMode:  "<<commonMode<<endl;
    /*
    c->cd(1);
    hh->Draw();
    c->cd(2);
    h->Draw();
    */

    // PERFORM COMMON MODE CORRECTION FOR A GIVEN TIME BIN
    std::transform (rawDataTS_1stSet.begin(), rawDataTS_1stSet.end(), rawDataTS_1stSet.begin(), std::bind2nd(std::minus<Float_t>(), commonMode));

    //  ADC SUM OVER ALL TIME BINS USE AS THE TEST CRITERIA FOR ZERO SUPPRESSION
    if(timebin == 0)  rawDataZS_1stSet.resize(rawDataTS_1stSet.size()) ;
    std::transform (rawDataZS_1stSet.begin(), rawDataZS_1stSet.end(), rawDataTS_1stSet.begin(), rawDataZS_1stSet.begin(), std::plus<Float_t>());
 
    //======================================================================
    //                    2ND SET OF CHANNELS
    //======================================================================

    // PERFORM APV25 PEDESTAL OFFSET CORRECTION  FOR A GIVEN TIME BIN
    std::transform (rawDataTS_2ndSet.begin(), rawDataTS_2ndSet.end(), fPedestalOffsets_2ndSet.begin(), rawDataTS_2ndSet.begin(), std::minus<Float_t>());

    // COMPUTE THRESOLD (= 20 PEDESTAL NOISE) FOR COMMON MODE DATA SELECTION 
    threshold_2ndSet.clear();
    threshold_2ndSet.resize( rawDataTS_2ndSet.size()  );
    std::transform (fPedestalNoises_2ndSet.begin(), fPedestalNoises_2ndSet.end(), threshold_2ndSet.begin(), std::bind1st(std::multiplies<Float_t>(), 20.0)); 

    // TEST IF APV CHANNEL IS A HIT OR PEDESTAL => ROUGH TEST TO ELIMINATE HIT CHANNELS IN FOR THE COMMON MODE
    dataTest_2ndSet.clear();
    dataTest_2ndSet.insert(dataTest_2ndSet.begin(), rawDataTS_2ndSet.begin(), rawDataTS_2ndSet.end());

    int nbHits_2ndSet = 0;
    for(int i=0;i<size2;i++)
    {
      if( (dataTest_2ndSet[i] - threshold_2ndSet[i]) >= 0)
      {
        dataTest_2ndSet[i] = 0;
	nbHits_2ndSet++;
      }
    }
    assert( nbHits_2ndSet <= size2);

    // COMPUTE COMMON MODE AND PERFORM COMMON MODE CORRECTION FOR A GIVEN APV AND TIME BIN FOR 2ND SET OF CHANNELS
    commonMode = std::accumulate(dataTest_2ndSet.begin(),dataTest_2ndSet.end(), 0.0) / (112 - nbHits_2ndSet) ;
    commonModeOffsets_2ndSet.push_back(commonMode) ;
    std::transform (rawDataTS_2ndSet.begin(), rawDataTS_2ndSet.end(), rawDataTS_2ndSet.begin(), std::bind2nd(std::minus<Float_t>(), commonMode));

    //  ADC SUM OVER ALL TIME BINS USE AS THE TEST CRITERIA FOR ZERO SUPPRESSION
    if(timebin == 0){ rawDataZS_2ndSet.clear(); rawDataZS_2ndSet.resize( rawDataTS_2ndSet.size()  );  }
    std::transform (rawDataZS_2ndSet.begin(), rawDataZS_2ndSet.end(), rawDataTS_2ndSet.begin(), rawDataZS_2ndSet.begin(), std::plus<Float_t>());
    
    //===============================================================================
    //        CLEAR EVERYTHING
    //===============================================================================
    rawDataTS.clear() ;
    rawDataTS_1stSet.clear() ;
    rawDataTS_2ndSet.clear() ;
    threshold_1stSet.clear() ;
    threshold_2ndSet.clear() ;
    dataTest_1stSet.clear() ;
    dataTest_2ndSet.clear() ;
 
    // PROCEED TO NEXT TIME BIN
    firstdata = lastdata + 12 ;
    lastdata = firstdata + NCH ;
    /*
    c->Update();
    getchar();
    c->Delete();
    h->Delete();
    */
  }
  //getchar();


  // ADC AVERAGE OVER ALL TIME BINS USE AS THE TEST CRITERIA FOR ZERO SUPPRESSION
  std::transform (rawDataZS_1stSet.begin(), rawDataZS_1stSet.end(), rawDataZS_1stSet.begin(), std::bind2nd(std::divides<Float_t>(), nbOfTimeBin)) ;
  std::transform (rawDataZS_2ndSet.begin(), rawDataZS_2ndSet.end(), rawDataZS_2ndSet.begin(), std::bind2nd(std::divides<Float_t>(), nbOfTimeBin)) ;

  //=============================================================================================================================
  // BLOCK WHERE ZERO SUPPRESSION IS COMPUTED
  //=============================================================================================================================
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
	avgdata = -rawDataZS_1stSet[i2nd] ;
	i2nd++;
      }

      // PERFORM ZERO SUPPRESSION
      if( fZeroSupCut > 0 ) {
	if( avgdata > (fZeroSupCut * fPedestalNoises[chNo]) ) {
	  if(!fListOfHits[hitID]) {
	    GEMHit * hit = new GEMHit(hitID, fAPVID, chNo, fZeroSupCut,  fIsHitMaxOrTotalADCs) ;
	    fListOfHits[hitID] = hit ;
	  }
	  fListOfHits[hitID]->AddTimeBinADCs(timebin, data) ;
	}
      }

      // NO ZERO SUPPRESSION
      else {
	Float_t adc = data ;
	if( avgdata < (5 * fPedestalNoises[chNo]) ) adc = avgdata ;
	if(!fListOfHits[hitID]) {
	  GEMHit * hit = new GEMHit(hitID, fAPVID, chNo, fZeroSupCut,  fIsHitMaxOrTotalADCs) ;
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
  threshold_1stSet.clear() ;
  threshold_2ndSet.clear() ;
  dataTest_1stSet.clear() ;
  dataTest_2ndSet.clear() ;

  rawDataTS.clear() ;
  commonModeOffsets_1stSet.clear() ;
  commonModeOffsets_2ndSet.clear() ;
}


