#ifndef GEM_EVENT_ANALYZER_H
#define GEM_EVENT_ANALYZER_H

#include <unordered_map>
#include <vector>
#include <string>

class GEMRawDecoder;
class GEMMapping;

class GEMEventAnalyzer
{
public:
    GEMEventAnalyzer();
    ~GEMEventAnalyzer();
    void Process(GEMRawDecoder * gem_decoder);
    void Clear();
    void SetEventNumber(unsigned int &);
    unsigned int & GetEventNumber();
    void ShowFecEvent(GEMRawDecoder * gem_decoder);
    int newEvent();
    void SetRunType(std::string);
    GEMRawDecoder * GetGEMRawDecoder();
    std::unordered_map<int, std::vector<int> > & GetEvent();
    bool IsEventUpdated();

private:
    GEMRawDecoder* gem_decoder;
    GEMMapping *mapping;
    unsigned int eventNumber;
    unsigned int lastEvent;
    bool new_event;
    unsigned int nFec;
    std::unordered_map<int, std::vector<int> > gem_event;

    // run type
    int physics_run;
    int pedestal_run;
    int raw_run;
};

#endif
