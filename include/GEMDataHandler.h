#ifndef GEMDATAHANDLER_H
#define GEMDATAHANDLER_H

#include "GEMDataStruct.h"

class GEMEventAnalyzer;
class GEMRawDecoder;
class EpicsEventAnalyzer;
class TDCEventAnalyzer;

class GEMDataHandler
{
public:
    GEMDataHandler();
    ~GEMDataHandler();
    void Decode(unsigned int *, int);
    void DecodeEpics(unsigned int *, int);
    void DispatchByROC(unsigned int *, int &);
    void DispatchByBank(unsigned int *, int &);
    GEMRawDecoder* GetGEMRawDecoder();
    void SetGEMEventAnalyzer(GEMEventAnalyzer * gem_event_analyzer);
    void SetEpicsEventAnalyzer(EpicsEventAnalyzer *fepics_analyzer);
    void SetTDCEventAnalyzer(TDCEventAnalyzer * ftdc_analyzer);

private:
    GEMRawDecoder *gem_decoder; 
    GEMEventAnalyzer *gemEventAnalyzer;
    EventHeader *headerROC;
    EventHeader *headerBank;
    EpicsEventAnalyzer * epics_analyzer;
    TDCEventAnalyzer * tdc_event_analyzer;

    unsigned int eventNumber;
};

#endif
