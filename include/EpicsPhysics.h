#ifndef EPICS_PHYSICS_H
#define EPICS_PHYSICS_H
#include <unordered_map>
#include <string>

class GEMTree;

class EpicsPhysics
{
public:
    EpicsPhysics();
    ~EpicsPhysics();

    void AccumulateEvent(unsigned int &, std::unordered_map<std::string, double> &);
    void SaveEpicsResults();
    void SetGEMTree(GEMTree *);

private:
    GEMTree * epics_tree;
};

#endif
