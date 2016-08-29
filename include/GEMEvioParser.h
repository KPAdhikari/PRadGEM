#ifndef GEMEVIOPARSER_H
#define GEMEVIOPARSER_H
#include <string>

#define BUF_SIZE 1000000

class GEMDataHandler;
class GEMConfigure;
class PRadDataHandler;
class EventUpdater;

class GEMEvioParser
{
public:
    GEMEvioParser();
    ~GEMEvioParser();

    void ParseFile(std::string &);
    void ParseEvent(unsigned int *);
    void SetDataHandler( GEMDataHandler * fHandler);
    void SetPRadDataHandler( PRadDataHandler *pHandler);
    void SetEventUpdater(EventUpdater *);
private:
    GEMDataHandler *handler;
    GEMConfigure * configure;
    PRadDataHandler *pHandler;
    EventUpdater *update_event;
    int eventLimit;
    int limit;
};

#endif
