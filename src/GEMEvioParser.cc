#include "GEMEvioParser.h"
#include "GEMDataHandler.h"
#include "GEMDataStruct.h"
#include "datastruct.h"
#include "GEMConfigure.h"
#include <evioUtil.hxx>
#include <evioFileChannel.hxx>
#include "PRadDataHandler.h"
#include "EventUpdater.h"

#define HEADER_SIZE 2

using namespace std;
using namespace evio;

GEMEvioParser::GEMEvioParser() {
    eventLimit = 3000;
    limit = 0;
    configure = new GEMConfigure();
    configure -> LoadConfigure();
    eventLimit = configure->GetNumEvtForPed();
}

GEMEvioParser::~GEMEvioParser()
{
}

void GEMEvioParser::SetDataHandler(GEMDataHandler *fHandler)
{
    handler = fHandler;
}

void GEMEvioParser::SetEventUpdater(EventUpdater *fupdate)
{
    update_event = fupdate;
}

void GEMEvioParser::SetPRadDataHandler(PRadDataHandler * handler)
{
    pHandler = handler;
}

void GEMEvioParser::ParseFile(string & file)
{
    cout<<"evio Parser: pasing file: "
	<<file<<endl;
    try {
	if(limit == eventLimit)
	    return;
	evioFileChannel chan(file.c_str(),"r", BUF_SIZE);
	chan.open();
	int event_index = 0;

	while(chan.read())
	{
	    if(event_index == eventLimit)
		break;
	    event_index ++;
	    limit ++;
	    //cout<<"event: "<<event_index<<endl;
	    ParseEvent(const_cast<unsigned int *> (chan.getBuffer() ) );
	    // chao's PRadDataHandler
	    pHandler->Decode(chan.getBuffer());
	    update_event->Update();
	}
	chan.close();
    } catch (evioException e) {
	cerr << endl <<e.toString()<<endl<<endl;
    }

}

void GEMEvioParser::ParseEvent(unsigned int * buf)
{
    EventHeader * header = (EventHeader*) &buf[0];

    switch(header->tag)
    {
	case CODA_Event:
	    {
		handler -> Decode(&buf[HEADER_SIZE], buf[0]-HEADER_SIZE);
		break;
	    }
	case CODA_Sync:
	    break;
	case EPICS_Info:
	    {
		handler -> DecodeEpics(&buf[HEADER_SIZE], buf[0]-HEADER_SIZE);
		break;
	    }
	case CODA_Prestart:
	    break;
	case CODA_Go:
	    break;
	case CODA_End:
	    break;
	default:
	    return;
    }
}

