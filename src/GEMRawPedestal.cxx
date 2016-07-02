#include "GEMRawPedestal.h"

static string trim(const string &str, const string &w = " \t\n\r")
{

    const auto strBegin = str.find_first_not_of(w);
    if (strBegin == string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(w);

    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


GEMRawPedestal::GEMRawPedestal( map<int, map<int, vector<int> > > map)
{
  mapping = PRDMapping::GetInstance();
  mSrsSingleEvent.clear();
  mSrsSingleEvent = map;
  vSingleApvData.clear();
  mStripOffset.clear();
  mStripRMS.clear();

  fApvHeaderLevel = 1500.0;
  fTimeSample = 3;
  //fTimeSample = 12;
  fAPVStatus = "normal";
  fAPVID = 0xffff;
  NCH = 128;

  fFECID = 999;
  fADCCh = 999;

  ComputeEventPedestal();
}

GEMRawPedestal::~GEMRawPedestal()
{
  ClearMap();
  vSingleApvData.clear();
}

void GEMRawPedestal::ClearMap()
{
  map<int, map<int, vector<int> > >::iterator it;
  for(it = mSrsSingleEvent.begin(); it!=mSrsSingleEvent.end();++it)
  {
    map<int, vector<int> > fecdata = it->second;

    map<int, vector<int> >::iterator itt;
    for(itt = fecdata.begin(); itt!=fecdata.end(); ++itt)
    {
      itt->second.clear();
    }

    it->second.clear();
  }
  mSrsSingleEvent.clear();

  map<int, map<int, vector<Float_t> > >::iterator it1, IT1;
  for(it1 = mStripOffset.begin(); it1!=mStripOffset.end();++it1)
  {
    map<int, vector<Float_t> > fecdata = it1->second;
    map<int, vector<Float_t> >::iterator it2;
    for(it2 = fecdata.begin(); it2!=fecdata.end(); ++it2)
    {
      it2->second.clear();
    }
    it1->second.clear();
  }
  mStripOffset.clear();
  for(IT1 = mStripRMS.begin(); IT1!=mStripRMS.end();++IT1)
  {
    map<int, vector<Float_t> > fecdata = IT1->second;
    map<int, vector<Float_t> >::iterator IT2;
    for(IT2=fecdata.begin(); IT2!=fecdata.end(); ++IT2)
    {
      IT2->second.clear();
    }
    IT1->second.clear();
  }
  mStripRMS.clear();
}

void GEMRawPedestal::ComputeEventPedestal()
{
  map<int, map<int, vector<int> > >::iterator it;
  for(it=mSrsSingleEvent.begin(); it!=mSrsSingleEvent.end();++it)
  {
    int fecid = it->first;
    map<int, vector<int> > fecdata = it->second;
    map<int, vector<int> >::iterator itt;
    for(itt=fecdata.begin(); itt!=fecdata.end(); ++itt)
    {
      int adcch= itt->first;
      vector<int> apvdata = itt->second;
      //cout<<"-----------------------------------------------------"<<endl;
      //cout<<"FEC:  "<<fecid<<"   ADC Ch:  "<<adcch<<endl;
      vSingleApvData.clear();
      vSingleApvData = apvdata;
      int apvid = fecid << 4 | adcch;
      fAPVID = apvid;
      if( mapping->GetAPVHeaderLevelFromID(apvid) ) fApvHeaderLevel = mapping->GetAPVHeaderLevelFromID(apvid); 
      
      fAPVStatus = mapping->GetAPVstatus(apvid); 
      fFECID = fecid;
      fADCCh = adcch;
      ApvEventDecode();
      
      ComputeApvPedestal(fecid, adcch);
      //cout<<"_____________________________________________________"<<endl;
    }
  }
}

void GEMRawPedestal::ApvEventDecode()
{
  Int_t idata = 0;
  Int_t startDataFlag = 0;

  int size = vSingleApvData.size();
  vCommonModeOffset.clear();
  vCommonModeOffset_split.clear();
  mApvTimeBinData.clear();

  std::list<Float_t> commonModeOffset;
  std::list<Float_t> commonModeOffset_split;

  //debug
  int nTB = 0;
  vector<int> iheader;
  iheader.clear();

  //----------------------------------------------------------//
  //  check time samples, if not right, discard current data  //
  //----------------------------------------------------------//
  int i = 0;
  int nbTimeSample = 0;
  while( i < size )
  {
    //cout<<i<<endl;
    if( vSingleApvData[i] < fApvHeaderLevel )
    {
      i++;
      if( vSingleApvData[i] < fApvHeaderLevel )
      {
        i++;
        if( vSingleApvData[i] < fApvHeaderLevel )
        {
          if(i+138<size)
          {
            i+=10; //8 address + 1 error bit
            i+=128;
            nbTimeSample++;
            iheader.push_back(idata);
            continue;
          }
        }
      }
    }
    i++;
  }
  //cout<<nbTimeSample<<endl;
  if( nbTimeSample != fTimeSample )
  {
    cout<<"GEMRawPedestal ERROR: Time Bins expected: "<<fTimeSample<<"  Time Bins computed: "<<nbTimeSample<<endl;
    cout<<"please"<<endl;
    cout<<"      1, check TS value in configure file."<<endl;
    cout<<"      2, or adjust APV Header in your mapping file."<<endl;
    cout<<"      3, or check raw data, make sure all APVs are working."<<endl;
    cout<<"      4, size: "<<size<<"  idata: "<<idata<<endl;
    cout<<"      5, FECID: "<<fFECID<<"  ADC_CH:  "<<fADCCh<<endl;
    cout<<"      6, Hit Enter to ignore this data sample ..."<<endl;
    int s = iheader.size();
    for(int ii = 0;ii<s;ii++)
    {
      cout<<iheader[ii]<<"  ";
    }
    cout<<endl;

    TCanvas *c = new TCanvas("c", "c", 10, 10, 600,400);
    TH1F* h = new TH1F("h", Form("FECID: %d, ADC_CH: %d", fFECID, fADCCh), size, 0, size);
    for(int i=0;i<size;i++)
    {
      h->SetBinContent(i+1, vSingleApvData[i]);
    }
    h->SetStats(0);
    h->Draw();
    c->Update();
    getchar();
    h->Delete();
    c->Delete();
  }

  if (nbTimeSample != fTimeSample) 
  {
    return;
  }
  //----------------------------------------------
  //  end check, ignore incorrect apvs...
  //----------------------------------------------

  while( idata < size )
  {
    //---------------------- --------------------------------//
    //   look for apv header  => 3 consecutive words < header//
    //   ----------------------------------------------------//
    if( vSingleApvData[idata] < fApvHeaderLevel )
    {
      idata++;
      if( vSingleApvData[idata] < fApvHeaderLevel )
      {
        idata++;
	if( vSingleApvData[idata] < fApvHeaderLevel )
	{
	  if(idata+138<size)
	  {
	    idata+=10; //8 address + 1 error bit
	    startDataFlag = 1;
	    nTB++;
	    iheader.push_back(idata);
	    continue;
	  }
	}
      }
    }

    //--------------------------------------
    //       meaningful data
    //--------------------------------------
    if(startDataFlag)
    {
      commonModeOffset.clear();
      commonModeOffset_split.clear();
      Float_t commonMode = 0;
      Float_t commonMode_split = 0;

      for(int chNo = 0; chNo < NCH; ++chNo)
      {
        Int_t stripNo = chNo;
        Float_t rawdata = vSingleApvData[idata];
	string apv_status = fAPVStatus.Data();
	apv_status = trim(apv_status);
	if( (apv_status != "normal") && (mapping->GetPRadStripMapping(fAPVID, chNo)<16) )
	{
          //cout<<"apv event decoder: not normal apvs:  "<<apv_status<<endl;
	  commonModeOffset_split.push_back(rawdata);
	}
	else
	{
          //cout<<"apv event decoder: normal apvs:  "<<apv_status<<endl;
	  commonModeOffset.push_back(rawdata);
	}
	mApvTimeBinData.insert(pair<Int_t, Float_t>(stripNo, rawdata) );
	idata++;
      }

      commonModeOffset.sort();
      commonModeOffset_split.sort();
      commonMode = TMath::Mean(commonModeOffset.begin(), commonModeOffset.end() );
      commonMode_split = TMath::Mean(commonModeOffset_split.begin(), commonModeOffset_split.end());
      commonModeOffset.clear();
      commonModeOffset_split.clear();

      vCommonModeOffset.push_back(commonMode);
      vCommonModeOffset_split.push_back(commonMode_split);
      startDataFlag = 0;
      continue;
    }
    idata++;
  }
  //debug
  //cout<<"APVEventDecode, Time bins found:  "<<nTB<<endl;
 
}

void GEMRawPedestal::ComputeApvPedestal(int fecid, int adcch)
{
  // compute offset and rms for each strip
  vector<Float_t> timeBinPedestalOffset;
  vector<Float_t> timeBinPedestalNoise;

  pair<multimap<Int_t, Float_t>::iterator, multimap<Int_t, Float_t>::iterator> stripTimeBinRawData;
  for(int stripNo = 0; stripNo<NCH; stripNo++)
  {   

    stripTimeBinRawData = mApvTimeBinData.equal_range(stripNo);
   
    multimap<Int_t, Float_t>::iterator timebin_it;
    Int_t timebin = 0;
    for(timebin_it = stripTimeBinRawData.first; timebin_it != stripTimeBinRawData.second; ++timebin_it)
    {
      Float_t rawdata = timebin_it->second;
      Float_t commonModeOffset = vCommonModeOffset[timebin];

      string apv_status = fAPVStatus.Data();
      apv_status = trim(apv_status);

      if( (apv_status != "normal") && (mapping->GetPRadStripMapping(fAPVID,stripNo)<16) )
      {
        //cout<<" compute_pedestal function:  not normal apvs:  "<<apv_status<<endl;
        commonModeOffset = vCommonModeOffset_split[timebin];
      }
      //timeBinPedestalNoise.push_back(rawdata - vCommonModeOffset[timebin] );
      timeBinPedestalNoise.push_back(rawdata - commonModeOffset );
      timeBinPedestalOffset.push_back(rawdata);

      timebin++;
    }

    Float_t meanPedestalNoise = TMath::Mean(timeBinPedestalNoise.begin(), timeBinPedestalNoise.end() );
    Float_t meanPedestalOffset = TMath::Mean(timeBinPedestalOffset.begin(), timeBinPedestalOffset.end() );

    // Kondo Way
    mStripOffset[fecid][adcch].push_back(meanPedestalOffset);
    mStripRMS[fecid][adcch].push_back( meanPedestalNoise);

    timeBinPedestalNoise.clear();
    timeBinPedestalOffset.clear();
  } 

}

Float_t GEMRawPedestal::GetStripNoise(int fecid, int adc_ch, int channelID)
{
  vector<Float_t> fvec = mStripRMS[fecid][adc_ch];
  assert(fvec.size() == NCH );
  return fvec[channelID];
}

Float_t GEMRawPedestal::GetStripOffset(int fecid, int adc_ch, int channelID)
{
  vector<Float_t> fvec = mStripOffset[fecid][adc_ch];
  assert(fvec.size() == NCH );
  return fvec[channelID];
}

// offset after common mode subtraction, useless.
Float_t GEMRawPedestal::GetStripMeanOffset(int fecid, int adc_ch, int channelID)
{
  vector<Float_t> fvec = mStripRMS[fecid][adc_ch];
  assert(fvec.size() == NCH );
  return fvec[channelID];
}

//debug member function
void GEMRawPedestal::PrintEventPedestal()
{
  map<int, map<int, vector<int> > >::iterator it;
  for(it = mSrsSingleEvent.begin(); it!=mSrsSingleEvent.end(); ++it)
  {
    int fecid = it->first;
    map<int, vector<int> > fecdata = it->second;

    map<int, vector<int> >::iterator itt;
    for(itt = fecdata.begin(); itt!=fecdata.end(); ++itt)
    {
      int adcch = itt->first;
      cout<<" ======================================"<<endl;
      cout<<"FEC:  "<<fecid<<"  adc_ch:  "<<adcch<<endl;
      for(int i=0;i<NCH;i++)
      {
        cout<<i<<":  Noise: "<<GetStripNoise(fecid, adcch, i);
	cout<<":  Offset:"<<GetStripOffset(fecid, adcch, i)<<endl;
      }
      cout<<" ======================================"<<endl;
    }
  }
}
