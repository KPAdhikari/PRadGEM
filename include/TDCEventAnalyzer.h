#ifndef TDC_EVENT_ANALYZER_H
#define TDC_EVENT_ANALYZER_H

#include <unordered_map>
#include <string>
#include <vector>

class TDCEventAnalyzer
{
public: 
    TDCEventAnalyzer();
    ~TDCEventAnalyzer();
    void ProcessEvent(unsigned int *, int);
    void PrintBuffer();
    void InitTDCGroup();
    void Decode();
    void SetEventNumber(unsigned int );
    unsigned int GetEventNumber();
    std::vector<double> GetTDCChannelValue(std::string);
    std::unordered_multimap<std::string, double> & GetEvent();

private:
    std::unordered_map<std::string, unsigned short> tdc_map;
    std::unordered_map<unsigned short, std::string> inverse_tdc_map;
    std::unordered_multimap<std::string, double> channel_value_map;
    // use multimap because for each tdc channel, 
    // there probably be multihits
    unsigned int * buf;
    int size;
    unsigned int event_number;
};

#endif
