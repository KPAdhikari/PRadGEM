#ifndef GEMANALYZER_H
#define GEMANALYZER_H

#include <string>

class GEMConfigure;
class GEMEvioParser;
class GEMDataHandler;
class GEMEventAnalyzer;
class GEMPedestal;
class GEMPhysics;
class PRadDataHandler;
class PRadReconstructor;
class EpicsEventAnalyzer;
class EpicsPhysics;
class GEMTree;
class TDCEventAnalyzer;
class GEMEfficiency;
class EventUpdater;

class GEMAnalyzer
{
public:
    GEMAnalyzer();
    ~GEMAnalyzer();

    void ProcessFiles();
    void Process();
    void SetDataHandler();
    void RegisterClass();
    void ProcessResults();
    void GetRunType();

    void InitPRadDataHandler();

private:
    std::string file;
    GEMConfigure *configure;
    GEMEvioParser *parser;
    GEMDataHandler *handler;
    GEMEventAnalyzer *gem_event_analyzer;
    GEMPedestal * pedestal;
    GEMPhysics * physics;
    EpicsEventAnalyzer * epics_event_analyzer;
    EpicsPhysics * epics_physics;
    GEMTree *res_tree;
    TDCEventAnalyzer *tdc_event_analyzer;
    GEMEfficiency * gem_efficiency;
    EventUpdater * event_updater;

    // chao
    PRadDataHandler * pHandler;
    PRadReconstructor *reconstruct;

    // run type
    int gem_physics;
    int gem_pedestal;
    int gem_raw;
};

#endif
