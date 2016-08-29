#include "GEMDataHandler.h"
#include "GEMRawDecoder.h"
#include "GEMDataStruct.h"
#include "datastruct.h"
#include "GEMEventAnalyzer.h"
#include "EpicsEventAnalyzer.h"
#include "TDCEventAnalyzer.h"
#include <iostream>
#include <iomanip>

#define HEADER_SIZE 2

using namespace std;

GEMDataHandler::GEMDataHandler()
{
    gem_decoder = new GEMRawDecoder();
}

GEMDataHandler::~GEMDataHandler()
{
}

void GEMDataHandler::SetGEMEventAnalyzer(GEMEventAnalyzer* gem_event_analyzer)
{
    gemEventAnalyzer = gem_event_analyzer;
}

void GEMDataHandler::Decode(unsigned int * buf, int size)
{
    int index = 0;

    // top level event id bank
    if( buf[index] != 4 )
    {
	cout<<"GEMDataHandler: Missing event id bank..."
	    <<endl;
	return;
    }
    eventNumber = buf[index+2];
    index += buf[index] +1;

    while(index < size)
    {
	// top level ROC bank
	headerROC = (EventHeader *) & buf[index];
	int len = buf[index] - HEADER_SIZE;
	DispatchByROC( &buf[index+HEADER_SIZE], len);
	index += buf[index] +1;
    }
}

void GEMDataHandler::DecodeEpics(unsigned int * buf, int size)
{
    int index = 0;
    while(index < size)
    {
        headerROC = (EventHeader *)&buf[index];
	int len = buf[index] - HEADER_SIZE;
	DispatchByROC( &buf[index+HEADER_SIZE], len);
	index += buf[index]+1;
    }
}

void GEMDataHandler::DispatchByROC(unsigned int *buf, int &size)
{
    int index = 0;
    while(index < size)
    {
	headerBank = (EventHeader *) & buf[index];
	int len = buf[index] -1;
	DispatchByBank(&buf[index + HEADER_SIZE], len);
	index += buf[index] +1;
    }
}

void GEMDataHandler::DispatchByBank(unsigned int *buf, int & size)
{
    switch(headerBank->tag)
    {
	case TI_BANK:
            break; 
	case EPICS_BANK:
        {
	    epics_analyzer -> SetEventNumber(eventNumber);
	    epics_analyzer -> ReadBuffer(buf);
	    //epics_analyzer -> PrintEpics();
	    epics_analyzer -> ProcessEvent();
	    break;
	}
	case GEM_BANK:
        {
            //cout<<"bank: "<<headerBank->tag<<endl;
            //cout<<"GEM BANK Size: "<<size<<endl;
	    gem_decoder->DecodeFEC(buf, size);
	    gemEventAnalyzer->SetEventNumber(eventNumber);
            gemEventAnalyzer->Process(gem_decoder);
	    break;
	}
	case TDC_BANK:
	{
	    tdc_event_analyzer->SetEventNumber(eventNumber);
	    tdc_event_analyzer->ProcessEvent(buf, size);
	    break;
	}
	case EVINFO_BANK:
	break;
	default:
	return;
    }
}

GEMRawDecoder* GEMDataHandler::GetGEMRawDecoder()
{
    return gem_decoder;
}

void GEMDataHandler::SetEpicsEventAnalyzer(EpicsEventAnalyzer * fepics_analyzer)
{
    epics_analyzer = fepics_analyzer;
}

void GEMDataHandler::SetTDCEventAnalyzer(TDCEventAnalyzer *fanalyzer)
{
    tdc_event_analyzer = fanalyzer;
}
