#ifndef GEMPHYSICS_H
#define GEMPHYSICS_H

#include <unordered_map>
#include <vector>
#include "datastruct.h"
#include "PRadEventStruct.h"

class GEMPedestal;
class GEMMapping;
class GEMOnlineHitDecoder;
class GEMTree;
class GEMCoord;
class PRadReconstructor;
class PRadDataHandler;
class PRadMoller;
class PRadEP;
class HyCalGEMMatch;

class GEMPhysics
{
public:
    GEMPhysics();
    ~GEMPhysics();
    void AccumulateEvent(int evtID, std::unordered_map<int, std::vector<int> > & event);
    void CharactorizeGEM();
    void CharactorizePhysics();
    void SavePhysResults();
    void SetPRadReconstructor(PRadReconstructor *reconstructor);
    void SetPRadDataHandler(PRadDataHandler * handler);
    void SetGEMTree(GEMTree *);

private:
    GEMPedestal *pedestal;
    GEMMapping * mapping;
    GEMOnlineHitDecoder * hit_decoder;
    GEMTree *rst_tree;
    GEMCoord *gem_coord;
    PRadMoller * moller_analyzer;
    PRadMoller * gem1_moller_analyzer; // for moller on single gem
    PRadMoller * gem2_moller_analyzer; // for moller on single gem
    PRadEP * ep_analyzer;
    HyCalGEMMatch * match_analyzer;

    PRadDataHandler *pHandler;
    PRadReconstructor * reconstruct;
    std::vector<HyCalHit> hycal_hit;
};

#endif
