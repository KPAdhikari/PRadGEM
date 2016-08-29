Target  := main

Source  := main.cc \
           GEMAnalyzer.cc \
	   GEMDataHandler.cc \
	   GEMMapping.cc \
	   GEMConfigure.cc \
	   GEMEvioParser.cc \
	   GEMRawDecoder.cc \
           PRadBenchMark.cc \
	   GEMEventAnalyzer.cc \
           GEMHit.cc \
           GEMCluster.cc \
           GEMRawPedestal.cc \
           GEMPedestal.cc \
           GEMPhysics.cc \
           GEMOnlineHitDecoder.cc \
           GEMTree.cc \
           GEMCoord.cc \
           EpicsEventAnalyzer.cc \
	   EpicsPhysics.cc \
	   HyCalGEMMatch.cc \
	   PRadMoller.cc \
	   PRadEP.cc \
	   TDCEventAnalyzer.cc \
	   GEMEfficiency.cc \
	   EventUpdater.cc \
	   MollerGEMSpatialRes.cc

OBJDIR 	:= ./obj

Objs	:= $(patsubst %.cc, %.o, $(Source))
OBJS	:= $(patsubst %.o, $(OBJDIR)/%.o, $(Objs))

cc      := g++

libs    := $(shell root-config --libs)
glibs   := $(shell root-config --glibs)
cflags  := $(shell root-config --cflags)

incfile := -I/home/xbai/w/coda/common/include -I/home/xbai/w/pRad/source/PRadDecoder/lib/include -I./include

flags   := -O3 -std=c++11 $(glibs) $(cflags) $(incfile) -L/home/xbai/w/coda/Linux-x86_64/lib -levio -levioxx -lexpat -L/home/xbai/w/pRad/source/PRadDecoder/lib -lPRadDecoder

$(Target) : $(OBJS)
	@$(cc) -o $(Target) $(OBJS) $(flags)

$(OBJDIR)/%.o: ./src/%.cc
	@echo Compiling $< ...
	@$(cc) -c $< -o $@ $(flags)

clean:
	@rm -f $(Target)
	@rm -f $(OBJDIR)/*.o
