//==============================================================//
// Process GEM Data                                             //
// Check GEM Raw                                                //
//                                                              //
// Xinzhan Bai                                                  //
// xb4zp@virginia.edu                                           //
//==============================================================//

#include "GEMEventAnalyzer.h"
#include "GEMRawDecoder.h"
#include "GEMMapping.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <iostream>

using namespace std;

GEMEventAnalyzer::GEMEventAnalyzer()
{
    eventNumber = 0;
    lastEvent = 0;
    mapping = GEMMapping::GetInstance();
    gem_event.clear();
    nFec = 0;
}

GEMEventAnalyzer::~GEMEventAnalyzer()
{
}

void GEMEventAnalyzer::SetRunType( string runtype)
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
        cout<<" GEMEventAnalyzer: undefined run type..."
	    <<endl;
	return;
    }
}

void GEMEventAnalyzer::Process(GEMRawDecoder * gem_decoder)
{
    if(raw_run)
    {
        ShowFecEvent(gem_decoder);
        return;
    }

    unordered_map<int, vector<int> > fecEvent = gem_decoder->GetFECDecoded();

    if(newEvent())
    {  
        gem_event.clear();
        lastEvent = eventNumber;
	gem_event.insert(fecEvent.begin(), fecEvent.end() );
    }
    else
    {
        gem_event.insert(fecEvent.begin(), fecEvent.end() );
        lastEvent = eventNumber;
    }
}

void GEMEventAnalyzer::SetEventNumber(unsigned int & n)
{
    eventNumber = n;
}

unsigned int & GEMEventAnalyzer::GetEventNumber()
{
    return eventNumber;
}

int GEMEventAnalyzer::newEvent()
{
    if( eventNumber  != lastEvent) {
        new_event = true;
	return 1;
    }
    else 
	return 0;
}

void GEMEventAnalyzer::ShowFecEvent(GEMRawDecoder * gem_decoder)
{
    cout<<"event number from data: "<<eventNumber<<endl;
    unordered_map<int, TH1F*> histo = gem_decoder->GetFECAPVRawHisto();
    TCanvas *c = new TCanvas("c", "c", 10, 10, 1200, 800);
    c->Divide(3,3);
    int ii=1;
    for(auto &i: histo)
    {
	c->cd(ii);
	i.second->Draw();
	ii++;
    }
    c->Update();
    getchar();
}

GEMRawDecoder * GEMEventAnalyzer::GetGEMRawDecoder()
{
    return gem_decoder;
}

unordered_map<int, vector<int> > & GEMEventAnalyzer::GetEvent()
{
    return gem_event;
}

bool GEMEventAnalyzer::IsEventUpdated()
{
    newEvent();
    if(new_event)
	return true;
    else 
	return false;
}
