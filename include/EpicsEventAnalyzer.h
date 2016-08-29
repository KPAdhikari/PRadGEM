#ifndef EPICSEVENTANALYZER_H
#define EPICSEVENTANALYZER_H

#include <string>
#include <queue>
#include <unordered_map>
#include <vector>

struct EpicsData
{
    double val;
    std::string name;

    EpicsData()
    : val(0), name("0")
    {}

    EpicsData(const double &v, const std::string &n)
    : val(v), name(n)
    {}
};

class EpicsEventAnalyzer
{
public:
    EpicsEventAnalyzer();
    ~EpicsEventAnalyzer();
    
    void ReadBuffer(unsigned int *);
    void PrintEpics();
    void ProcessEvent();
    void Strip();
    void SetEventNumber(const unsigned int &);
    unsigned int & GetEventNumber();
    std::unordered_map<std::string, double> & GetEpicsMap();
    bool IsEventUpdated();
    void newEvent();

private:
    unsigned int event_number;
    unsigned int last_event;
    bool new_event;
    std::queue<std::string> elements;
    std::vector<EpicsData> epics_vector;
    std::unordered_map<std::string, double> epics_map;
};

#endif
