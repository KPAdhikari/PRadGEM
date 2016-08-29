#ifndef __GEMPEDESTAL_H__
#define __GEMPEDESTAL_H__
#include <fstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <set>
#include <string>
#include <TH1F.h>
#include <TFile.h>

class GEMRawDecoder;
class GEMMapping;
class GEMRawPedestal;
class GEMConfigure;

class GEMPedestal
{
public:
    GEMPedestal();
    ~GEMPedestal();

    void BookHistos();
    void SetFile();

    int AccumulateEvent(int evtID, std::unordered_map<int, std::vector<int> > & event);
    void ComputePedestal();
    void SavePedestalFile();
    void LoadPedestal();
    void Delete();
    std::vector<Float_t> GetAPVNoises(int);
    std::vector<Float_t> GetAPVOffsets(int);

private:
    std::set<int> FECs;
    int NCH;
    int nNbofAPVs;

    std::string filename;

    ifstream file;

    std::vector<TH1F*> vStripOffsetHistos;
    std::vector<TH1F*> vStripNoiseHistos;
    std::vector<TH1F*> vApvPedestalOffset;
    std::vector<TH1F*> vApvPedestalNoise;

    TH1F* hAllStripNoise;
    TH1F* hAllXStripNoise;
    TH1F* hAllYStripNoise;

    // not care for now, will implement later if interested
    //TH1F* hAllStripOffset;
    //TH1F* hAllXStripOffset;
    //TH1F* hAllYStripOffset;

    GEMMapping * mapping;
    GEMConfigure *config;
    GEMRawPedestal *eventPedestal;
};

#endif
