#ifndef EVENT_UPDATER_H
#define EVENT_UPDATER_H
#include <string>

class GEMEventAnalyzer;
class TDCEventAnalyzer;
class GEMEfficiency;
class GEMPhysics;
class GEMPedestal;
class EpicsEventAnalyzer;
class EpicsPhysics;
class GEMTree;

class EventUpdater
{
public:
    EventUpdater();
    ~EventUpdater();

    void Update();
    //gem
    void ProcessGEMEfficiency();
    void ProcessGEM();
    void ProcessGEMPhysics();
    void ProcessGEMPedestal();
    void ProcessGEMRaw();
    void SetRunType(std::string runtype);
    //epics
    void ProcessEpics();
    void ProcessEpicsPhysics();

    void SetGEMEfficiencyAnalyzer(GEMEfficiency *);
    void SetGEMEventAnalyzer(GEMEventAnalyzer *);
    void SetTDCEventAnalyzer(TDCEventAnalyzer *);
    void SetGEMPhysics(GEMPhysics*);
    void SetGEMPedestal(GEMPedestal*);
    void SetGEMTree(GEMTree*);
    void SetEpicsEventAnalyzer(EpicsEventAnalyzer *);
    void SetEpicsPhysics(EpicsPhysics*);

private:
    GEMEfficiency * effi_analyzer;
    TDCEventAnalyzer *tdc_event_analyzer;
    GEMEventAnalyzer *gem_event_analyzer;
    GEMPhysics * physics;
    GEMPedestal* pedestal;
    GEMTree *tree;
    EpicsEventAnalyzer * epics_event_analyzer;
    EpicsPhysics * epics_physics;

    // gem run type
    int physics_run;
    int pedestal_run;
    int raw_run;
};

#endif
