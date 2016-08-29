#ifndef GEM_EFFICIENCY_H
#define GEM_EFFICIENCY_H
#include <vector>
#include <unordered_map>
#include <string>
#include "datastruct.h"
#include "PRadEventStruct.h"
#include <fstream>

#define TRACK_PERIOD 2000
#define CHECK_PAIR_COMPTON

class GEMConfigure;
class GEMTree;
class GEMCoord;
class GEMOnlineHitDecoder;
class GEMPedestal;
class GEMMapping;
class PRadReconstructor;
class PRadDataHandler;
class HyCalGEMMatch;

class GEMEfficiency
{
public:
    GEMEfficiency();
    ~GEMEfficiency();
    bool TDCCut(std::unordered_multimap<std::string, double> & tdc_map);
    bool HyCalTimingCut(std::unordered_multimap<std::string, double> & tdc_map);
    bool ScinTimingCut(std::unordered_multimap<std::string, double> & tdc_map);
    void InitTDCCutVar();

    void AccumulateEvent(unsigned int, std::unordered_map<int, std::vector<int> > &,
	    std::unordered_multimap<std::string, double> &);
    void IncreasePeriod(int);
    void IncreaseOverall(int);

    void SetGEMTree(GEMTree * tree);
    void SetEventNumber(unsigned int );
    unsigned int GetEventNumber();
    void SetGEMConfigure(GEMConfigure *);
    void SetGEMCoord(GEMCoord *);
    void SetPRadReconstructor(PRadReconstructor *);
    void SetPRadDataHandler(PRadDataHandler *);
    bool passHyCalGEMPosMatch();
    void SaveGEM();
    void InitLogFile();

    // align GEM 
    void AlignGEM();

private:
    GEMConfigure * configure;
    GEMTree * gem_tree;
    GEMCoord * gem_coord;
    GEMOnlineHitDecoder * hit_decoder;
    GEMPedestal * pedestal;
    GEMMapping * mapping;
    unsigned int event_number;

    PRadDataHandler * pHandler;
    PRadReconstructor * reconstructor;
    std::vector<HyCalHit> hycal_hit;

    HyCalGEMMatch * match_analyzer;

    // tdc cut 
    int hycal_group_q;
    std::string hycal_group[9];
    float hycal_start;
    float hycal_end;
    std::string scin_mode;
    int scin_q;
    std::string scin_group[2];
    float scin_start;
    float scin_end;

    // for eff. stability, every 2k events 
    // put down a eff. number
    std::vector<float> eff_track;
    float period_trigger;
    float period_gem;
    float n_trigger;
    float n_gem;
    std::fstream log_file;

#ifdef CHECK_PAIR_COMPTON
    // pair production, compton scattering event ratio check
    // single electron on gem1 && gem2
    double compton_event; 
    double total_event;
    // two electron event on gem1 || gem2
    double pair_event; 
#endif
};

#endif
