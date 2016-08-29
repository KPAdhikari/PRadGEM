#include "EpicsEventAnalyzer.h"
#include <iostream>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include "EpicsPhysics.h"

using namespace std;

EpicsEventAnalyzer::EpicsEventAnalyzer()
{
    event_number = 0;
    last_event = 0;
}

EpicsEventAnalyzer::~EpicsEventAnalyzer()
{
    queue<string>().swap(elements);
}

void EpicsEventAnalyzer::ReadBuffer(unsigned int *buf)
{
    queue<string>().swap(elements);
    string buffer = (char*)buf;
    string line;
    for(auto c : buffer)
    {
        if( c!= '\n')
	{
	    line += c;
	}
	else
	{
	    elements.push(line);
	    line ="";
	}
    }
}

void EpicsEventAnalyzer::PrintEpics()
{
    cout<<endl;
    cout<<"EpicsEventAnalyzer::PrintEpics..."
        <<endl<<endl;
    cout<<event_number<<endl;
    Strip();
    for(auto &i: epics_map)
    {
        cout<<setfill(' ')<<setw(16)<<i.first<<"  "<<setfill(' ')<<setw(16)<<i.second<<endl;
    }
}

void EpicsEventAnalyzer::ProcessEvent()
{
    Strip();
}

void EpicsEventAnalyzer::SetEventNumber(const unsigned int & n)
{
    event_number = n;
}

unsigned int & EpicsEventAnalyzer::GetEventNumber()
{
    return event_number;
}

void EpicsEventAnalyzer::Strip()
{
    epics_vector.clear();
    epics_map.clear();
    while(!elements.empty())
    {
	string s = elements.front();
	//cout<<s<<endl;
	elements.pop();
	if(s.empty())
	    continue;
	//s.erase( remove_if(s.begin(), s.end(), ::isspace), s.end());
	char *tokens = strtok( (char*)s.data(), " ");
	char *tok = strtok( nullptr, " ");
	epics_vector.emplace_back(atof(tokens), string(tok));
    }
    for(auto &i: epics_vector)
    {
	epics_map[i.name] = i.val;
    }
}

unordered_map<string, double> & EpicsEventAnalyzer::GetEpicsMap()
{
    return epics_map;
}

void EpicsEventAnalyzer::newEvent()
{
    if( event_number != last_event ) {
	new_event = true;
	last_event = event_number;
    }
    else 
	new_event = false;
}

bool EpicsEventAnalyzer::IsEventUpdated()
{
    newEvent();
    if( new_event )
	return true;
    else
	return false;

}
