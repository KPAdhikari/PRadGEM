#include "TDCEventAnalyzer.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

#define TDC_UNDEFINED -99999999

using namespace std;

TDCEventAnalyzer::TDCEventAnalyzer()
{
    InitTDCGroup();
}

TDCEventAnalyzer::~TDCEventAnalyzer()
{
}

void TDCEventAnalyzer::ProcessEvent(unsigned int * fbuf, int fsize)
{
    buf =  fbuf;
    size = fsize;
    Decode();
}

vector<double> TDCEventAnalyzer::GetTDCChannelValue(string str)
{
    vector<double> vec;
    auto range = channel_value_map.equal_range(str);
    if(range.first == range.second){
        vec.push_back(TDC_UNDEFINED);
	return vec;
    }

    for_each(
        range.first,
	range.second,
	[&vec](unordered_multimap<string, double>::value_type &x){
	    vec.push_back(x.second);
	}
    );

    return vec;
}

void TDCEventAnalyzer::Decode()
{
    channel_value_map.clear();
    for(int i=0;i<size;i++){
	if( (buf[i] & 0xf8000000) != 0)
	    continue;
	int tdc_ch = ( buf[i]>>19 ) & 0x7f; 
	if( inverse_tdc_map.find(tdc_ch) != inverse_tdc_map.end()) {
	    channel_value_map.insert(make_pair<string, double>( (string)inverse_tdc_map[tdc_ch], (double)(buf[i]&0x7ffff)));
	}
    }
}

void TDCEventAnalyzer::PrintBuffer()
{
    cout<<"TDC Print buffer."<<endl;
    cout<<"event number: "<<event_number<<endl;
    for(int i=0;i<size;i++)
    {
	cout<<setfill('0')<<setw(8)<<hex<<buf[i]<<endl;
    }
}

void TDCEventAnalyzer::SetEventNumber(unsigned int n)
{
    event_number = n;
}

unsigned int TDCEventAnalyzer::GetEventNumber()
{
    return event_number;
}

void TDCEventAnalyzer::InitTDCGroup()
{
    //--------------------------------------------
    // HyCal TDC Group 
    // G1  112  :   G2  100 :  G3  109 :  G4  64
    // G5  69   :   G6  113 :  G25 65  :  W6  66
    // G10 67   :   W11  68 :  W18 70  :  W15 72
    // W9  73   :   W20  74 :  W33 75  :  W8  76
    // W27 77   :   W21  78 :   W2 80  :  W26 81
    // W7  82   :   W1   83 :   W3 84  :  W14 85
    // W24 88   :   G24  89 :   G20 90 :  G15 91
    // W12 92   :   W36  93 :   W30 94 :  G11 96
    // G22 97   :   W25  98 :   W13 99 :  G2  100
    // W32 101  :   W28  104:   W4  105:  W10 106
    // W16 107  :   W22  108:   G3  109:  G23 110
    // G16 114  :  G21   115:   W31 116:  W19 117
    // W23 120  :   W34  121:   W17 122:  W35 123
    // W29 124  :    W5  125: Scin S1 126 : Scin S2 127
    //-------------------------------------------------
    tdc_map.clear();
    inverse_tdc_map.clear();
    tdc_map[string("G1")] = 112;tdc_map[string("G2")] = 100;
    tdc_map[string("G3")] = 109;tdc_map[string("G4")] = 64;
    tdc_map[string("G5")] = 69;tdc_map[string("G6")] = 113;
    tdc_map[string("G25")] = 65;tdc_map[string("W6")] = 66;
    tdc_map[string("G10")] = 67;tdc_map[string("W11")] = 68;
    tdc_map[string("W18")] = 70;tdc_map[string("W15")] = 72;
    tdc_map[string("W9")] = 73;tdc_map[string("W20")] = 74;
    tdc_map[string("W33")] = 75;tdc_map[string("W8")] = 76;
    tdc_map[string("W27")] = 77;tdc_map[string("W21")] = 78;
    tdc_map[string("W2")] = 80;tdc_map[string("W26")] = 81;
    tdc_map[string("W7")] = 82;tdc_map[string("W1")] = 83;
    tdc_map[string("W3")] = 84;tdc_map[string("W14")] = 85;
    tdc_map[string("W24")] = 88;tdc_map[string("G24")] = 89;
    tdc_map[string("G20")] = 90;tdc_map[string("G15")] = 91;
    tdc_map[string("W12")] = 92;tdc_map[string("W36")] = 93;
    tdc_map[string("W30")] = 94;tdc_map[string("G11")] = 96;
    tdc_map[string("G22")] = 97;tdc_map[string("W25")] = 98;
    tdc_map[string("W13")] = 99;tdc_map[string("G2")] = 100;
    tdc_map[string("W32")] = 101;tdc_map[string("W28")] = 104;
    tdc_map[string("W4")] = 105;tdc_map[string("W10")] = 106;
    tdc_map[string("W16")] = 107;tdc_map[string("W22")] = 108;
    tdc_map[string("G3")] = 109;tdc_map[string("G23")] = 110;
    tdc_map[string("G16")] = 114;tdc_map[string("G21")] = 115;
    tdc_map[string("W31")] = 116;tdc_map[string("W19")] = 117;
    tdc_map[string("W23")] = 120;tdc_map[string("W34")] = 121;
    tdc_map[string("W17")] = 122;tdc_map[string("W35")] = 123;
    tdc_map[string("W29")] = 124;tdc_map[string("W5")] = 125;
    tdc_map[string("S1")] = 126;tdc_map[string("S2")] = 127;

    for(auto &i:tdc_map)
	inverse_tdc_map[ i.second ] = i.first;
}

unordered_multimap<string, double> & TDCEventAnalyzer::GetEvent()
{
    return channel_value_map;
}
