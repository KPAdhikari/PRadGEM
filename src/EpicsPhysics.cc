#include "EpicsPhysics.h"
#include "GEMTree.h"
#include <iostream>

using namespace std;

EpicsPhysics::EpicsPhysics()
{
    cout<<"EpicsPhysics Constructor..."
        <<endl;
}

EpicsPhysics::~EpicsPhysics()
{
}

void EpicsPhysics::AccumulateEvent(unsigned int & evtID, unordered_map<string, double> & epics_map)
{
    epics_tree -> SetEventID(evtID);
    epics_tree -> PushEpics(epics_map);
    epics_tree -> FillEpicsTree();
}

void EpicsPhysics::SetGEMTree(GEMTree * tree)
{
    epics_tree = tree;
}

void EpicsPhysics::SaveEpicsResults()
{
    //epics_tree -> WriteToDisk();
}
