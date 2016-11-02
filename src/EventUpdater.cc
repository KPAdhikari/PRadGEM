#include "EventUpdater.h"
#include "GEMEfficiency.h"
#include "GEMEventAnalyzer.h"
#include "GEMTree.h"
#include "TDCEventAnalyzer.h"
#include "GEMPhysics.h"
#include "GEMPedestal.h"
#include "GEMRawDecoder.h"
#include <iostream>
#include <unordered_map>
#include <vector>
#include "EpicsEventAnalyzer.h"
#include "EpicsPhysics.h"

using namespace std;

#define EFFICIENCY_RUN

EventUpdater::EventUpdater()
{
}

EventUpdater::~EventUpdater()
{
}

void EventUpdater::Update()
{
    ProcessGEM();
    ProcessEpics();
    //ProcessTDC();
}

void EventUpdater::ProcessGEM()
{
#ifdef EFFICIENCY_RUN
    ProcessGEMEfficiency();
#else
    if(physics_run)
        ProcessGEMPhysics();
    else if(pedestal_run)
        ProcessGEMPedestal();
    else if(raw_run)
        return;
    else{
        cout<<"raw: "<<raw_run<<" physics: "<<physics_run<<" pedestal: "<<pedestal_run<<endl;
        cout<<"EventUpdater:: undefined gem run type..."
	    <<endl;
    } 
#endif
}

void EventUpdater::ProcessGEMEfficiency()
{
    if( ! gem_event_analyzer->IsEventUpdated() )
	return;

    int event_number = gem_event_analyzer->GetEventNumber();

    unordered_multimap<string, double> tdc_map = 
	tdc_event_analyzer->GetEvent();

    effi_analyzer->AccumulateEvent(event_number, gem_event_analyzer->GetEvent(), tdc_map);
}

void EventUpdater::ProcessGEMPhysics()
{
    if( ! gem_event_analyzer->IsEventUpdated() )
	return;
    int event_number = gem_event_analyzer->GetEventNumber();
    physics->AccumulateEvent(event_number, gem_event_analyzer->GetEvent());
}

void EventUpdater::ProcessGEMPedestal()
{
    if( ! gem_event_analyzer->IsEventUpdated() )
	return;
    int event_number = gem_event_analyzer->GetEventNumber();
    pedestal->AccumulateEvent(event_number, gem_event_analyzer->GetEvent());
}

void EventUpdater::ProcessEpics()
{
    ProcessEpicsPhysics();
}

void EventUpdater::ProcessEpicsPhysics()
{
    if( ! epics_event_analyzer->IsEventUpdated() )
	return;
    unsigned int event_number = epics_event_analyzer->GetEventNumber();
    //cout<<"epics event number: "<<event_number<<endl;
    epics_physics -> AccumulateEvent(event_number, epics_event_analyzer->GetEpicsMap());
}

void EventUpdater::SetGEMEfficiencyAnalyzer(GEMEfficiency * feff)
{
    effi_analyzer = feff;
}

void EventUpdater::SetGEMEventAnalyzer(GEMEventAnalyzer * analyzer)
{
    gem_event_analyzer = analyzer;
}

void EventUpdater::SetTDCEventAnalyzer(TDCEventAnalyzer * analyzer)
{
    tdc_event_analyzer = analyzer;
}

void EventUpdater::SetGEMTree(GEMTree * ftree)
{
    tree = ftree;
}

void EventUpdater::SetGEMPhysics(GEMPhysics * fphysics)
{
    physics = fphysics;
}

void EventUpdater::SetGEMPedestal(GEMPedestal * fpedestal)
{
    pedestal = fpedestal;
    if(pedestal_run)
	pedestal -> BookHistos();
}

void EventUpdater::SetEpicsEventAnalyzer(EpicsEventAnalyzer * analyzer)
{
    epics_event_analyzer = analyzer;
}

void EventUpdater::SetEpicsPhysics(EpicsPhysics * physics)
{
    epics_physics = physics;
}

void EventUpdater::SetRunType(string runtype)
{
    physics_run = 0;
    pedestal_run = 0;
    raw_run = 0;
    if(runtype == "PHYSICS")
	physics_run = 1;
    else if(runtype == "PEDESTAL")
	pedestal_run = 1;
    else if(runtype == "RAWDATA")
	raw_run = 1;
    else
    {
	cout<<"EventUpdater SetRunType: undefined run type..."
	    <<endl;
	return;
    }
}
