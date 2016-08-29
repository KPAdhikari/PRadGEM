#include "GEMRawPedestal.h"
#include "GEMMapping.h"
#include <TMath.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <iostream>
#include <cassert>

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

GEMRawPedestal::GEMRawPedestal()
{
    mapping = GEMMapping::GetInstance();
    mSrsSingleEvent.clear();
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
}

GEMRawPedestal::GEMRawPedestal(unordered_map<int, vector<int> > r_map)
{
    mapping = GEMMapping::GetInstance();
    mSrsSingleEvent.clear();
    mSrsSingleEvent = r_map;
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

void GEMRawPedestal::ComputeEventPedestal(unordered_map<int, vector<int> > r_map)
{
    ClearMap();
    vSingleApvData.clear();
    mSrsSingleEvent = r_map;
    ComputeEventPedestal();
}

GEMRawPedestal::~GEMRawPedestal()
{
    ClearMap();
    vSingleApvData.clear();
}

void GEMRawPedestal::ClearMap()
{
    ClearSingleEvent();
    ClearStripOffset();
    ClearStripRMS();
}

void GEMRawPedestal::ClearSingleEvent()
{
    for(auto &i: mSrsSingleEvent)
    {
        i.second.clear();
    }
    mSrsSingleEvent.clear();
}

void GEMRawPedestal::ClearStripOffset()
{
    for(auto &i: mStripOffset)
    {
        i.second.clear();
    }
    mStripOffset.clear();
}

void GEMRawPedestal::ClearStripRMS()
{
    for(auto &i : mStripRMS)
    {
        i.second.clear();
    }
    mStripRMS.clear();
}

void GEMRawPedestal::ComputeEventPedestal()
{
    for(auto &i : mSrsSingleEvent)
    {
	fAPVID = i.first;
	fADCCh = fAPVID & 0xf;
	fFECID = (fAPVID>>4)&0xf;
	fApvHeaderLevel = mapping->GetAPVHeaderLevelFromID(fAPVID);
	fAPVStatus = mapping->GetAPVstatus(fAPVID);
	vSingleApvData = i.second;
	ApvEventDecode();
	ComputeApvPedestal(fAPVID);
    }
}

void GEMRawPedestal::ApvEventDecode()
{
    if (! CheckApvTimeSample()) 
	return;

    int idata = 0;
    int startDataFlag = 0;
    int size = vSingleApvData.size();
    vCommonModeOffset.clear();
    vCommonModeOffset_split.clear();
    mApvTimeBinData.clear();
    string apv_status = fAPVStatus.Data();
    apv_status = trim(apv_status);
	
    std::list<float> commonModeOffset;
    std::list<float> commonModeOffset_split;

    int nTB = 0;
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
	    float commonMode = 0;
	    float commonMode_split = 0;

	    for(int chNo = 0; chNo < NCH; ++chNo)
	    {
		float rawdata = vSingleApvData[idata];
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
		mApvTimeBinData.insert(pair<int, float>(chNo, rawdata) );
		idata++;
	    }

	    commonModeOffset.sort();
	    commonModeOffset_split.sort();
	    commonMode = TMath::Mean(commonModeOffset.begin(), commonModeOffset.end() );
	    commonMode_split = TMath::Mean(commonModeOffset_split.begin(), commonModeOffset_split.end());

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

int GEMRawPedestal::CheckApvTimeSample()
{
    //----------------------------------------------------------//
    //  check time samples, if not right, discard current data  //
    //----------------------------------------------------------//
    vector<int> iheader;
    iheader.clear();
    int size = vSingleApvData.size();
 
    int i = 0;
    int idata = 0;
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
                        idata = i;
			iheader.push_back(i);
			i+=128;
			nbTimeSample++;
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
	cout<<"GEMRawPedestal ERROR: Time Bins expected: "
	    <<fTimeSample
	    <<"  Time Bins computed: "
	    <<nbTimeSample
	    <<endl;
	cout<<"please"
	    <<endl;
	cout<<"      1, check TS value in configure file."
	    <<endl;
	cout<<"      2, or adjust APV Header in your mapping file."
	    <<endl;
	cout<<"      3, or check raw data, make sure all APVs are working."
	    <<endl;
	cout<<"      4, size: "
	    <<size<<"  idata: "
	    <<idata
	    <<endl;
	cout<<"      5, FECID: "<<fFECID
	    <<"  ADC_CH:  "<<fADCCh
	    <<endl;
	cout<<"      6, Hit Enter to ignore this data sample ..."
	    <<endl;
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
	return 0;
    }
    return 1;
}

void GEMRawPedestal::ComputeApvPedestal(int apvindex)
{
    // compute offset and rms for each strip
    vector<float> timeBinPedestalOffset;
    vector<float> timeBinPedestalNoise;
    string apv_status = fAPVStatus.Data();
    apv_status = trim(apv_status);

    pair<multimap<int, float>::iterator, multimap<int, float>::iterator> stripTimeBinRawData;
    for(int stripNo = 0; stripNo<NCH; stripNo++)
    {   
	stripTimeBinRawData = mApvTimeBinData.equal_range(stripNo);

	multimap<int, float>::iterator timebin_it;
	int timebin = 0;
	for(timebin_it = stripTimeBinRawData.first; timebin_it != stripTimeBinRawData.second; ++timebin_it)
	{
	    float rawdata = timebin_it->second;
	    float commonModeOffset = vCommonModeOffset[timebin];

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

	float meanPedestalNoise = TMath::Mean(timeBinPedestalNoise.begin(), timeBinPedestalNoise.end() );
	float meanPedestalOffset = TMath::Mean(timeBinPedestalOffset.begin(), timeBinPedestalOffset.end() );

	// Kondo Way
	mStripOffset[apvindex].push_back(meanPedestalOffset);
	mStripRMS[apvindex].push_back( meanPedestalNoise);

	timeBinPedestalNoise.clear();
	timeBinPedestalOffset.clear();
    } 
}

float GEMRawPedestal::GetStripNoise(int apvindex, int channelID)
{
    assert(mStripRMS[apvindex].size() == NCH );
    return 
        mStripRMS[apvindex][channelID];
}

float GEMRawPedestal::GetStripOffset(int apvindex, int channelID)
{
    assert(mStripOffset[apvindex].size() == NCH );
    return mStripOffset[apvindex][channelID];
}

//debug member function
void GEMRawPedestal::PrintEventPedestal()
{
    for(auto &i: mSrsSingleEvent)
    {
        int apvid = i.first;
	int adcch = apvid & 0xf;
	int fecid = (apvid>>4) & 0xf;
	for(int i=0;i<NCH;i++)
	{
	    cout<<" ======================================"<<endl;
	    cout<<"FEC:  "<<fecid<<"  adc_ch:  "<<adcch<<endl;
	    for(int i=0;i<NCH;i++)
	    {
		cout<<i<<":  Noise: "<<GetStripNoise(apvid, i);
		cout<<":  Offset:"<<GetStripOffset(apvid, i)<<endl;
	    }
	    cout<<" ======================================"<<endl;
	}
    }
}
