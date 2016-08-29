
#section 1),
    VERSION1:
    Architechture:
---------------------------------------------------------------------------------------------------------------
    analysis class        |       detector class              |      data handler        |        evio parser
---------------------------------------------------------------------------------------------------------------       
pedestal analysis   ->|
 physics analysis   ->|register-> GEMEventAnalyer   ->register->|
                                                                |        
pedestal analysis   ->|                                         |
 physics analysis   ->|register-> HyCalEventAnalyer ->register->|-> GEMDataHandler  ->register-> GEMEvioParser
                                                                |
pedestal analysis   ->|                                         |
 physics analysis   ->|register-> EPICSEventAnalyer ->register->|
                                  
                     				  GEMAnalyzer
---------------------------------------------------------------------------------------------------------------
         All classes declared in GEMAnalyzer, and registration relationship happens in it too
---------------------------------------------------------------------------------------------------------------

    VERSION2:
    Architechture:
        a), EvioParser parse event to DataHandler;
        b), DataHandler dispatch bank data to different detector classes, e.g. GEMEventAnalyzer.
        c), EventUpdater get decoded data from different EventAnalyzer, and accumulate events.
        d), Pedestal or Physics be registered to EventUpdater, access events data here, 
	    (access EventAnalyzer already registered in EventUpdater). so after each event, if necessary,
            one detector's physics or pedestal class can have the chance to access other detector's
            pedestal or physics data.
	e), EventUpdater registered to EvioParser, so after each event, Anlysis results get updated.
    Registration Relationship:
        a), all registration happens in GEMAnalyzer class.
	b), DataHandler register to EvioParser.
	c), EventAnalyzer register to DataHandler.
	d), EventAnalyzer register to Physics and Pedestal....
        d'), EventAnalyzer register to EventUpdater.
	e), Pedestal and Physics register to EventUpdater.
	f), EventUpdater register to EvioParser.
    Dive into Detectors:
        A), GEM detectors:
	    -), RawDecoder register to EventAnalyzer
	        RawDecoder produce data from bank data
            -), OnlineHitDecoder register to EventAnalyzer
	        OnlineHitDecoder extracts hits from decoded data.

#section 2), 
File <datastruct.h>
    This file defines experimental setup, ROC IDs, Bank IDs
    are defined here, etc...

    For Kondo:
        In this file, FEC_BANK1 ~ 4 is the bank in use

File <GEMDataStruct.h>
    some basic data structure

File <GEMAnalyzer.h>
    This is an interface class, not doing real work,
    worker classes should be registered in here.
    (passed to GEMDataHandler class from here)

File <GEMEvioParser.h>
    This file parse evio events, it gets events,
    then pass it to GEMDataHandler

File< GEMDataHandler.h>
    This file Parse event data, it parse event by bank id,
    then according to different bank id(different detectors)
    , it calls the corresponding worker class, pass data to
    specific worker class

File <GEMEventAnalyzer.h>
    This file Process GEM data.
    worker class for GEM detector.

File <PRadBenchMark.h>
   borrowed...

File <GEMHit.h>
    you know this very well...
    No willing to bother...

File more to be added...

to be continued...

