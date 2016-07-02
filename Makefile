Target  := main

Source  := main.cxx \
           GEMRawDecoder.cxx \
	   GEMInputHandler.cxx \
	   PRDMapping.cxx \
	   GEMRawPedestal.cxx \
	   GEMPedestal.cxx \
	   GEMHit.cxx \
	   GEMHitDecoder.cxx \
	   GEMOnlineHitDecoder.cxx \
	   GEMCluster.cxx \
	   GEMZeroHitDecoder.cxx \
	   GEMConfigure.cxx \
	   GEMHistoManager.cxx \
	   GEMPhysHandler.cxx

OBJDIR 	:= ./obj

Objs	:= $(patsubst %.cxx, %.o, $(Source))
OBJS	:= $(patsubst %.o, $(OBJDIR)/%.o, $(Objs))

cc      := g++

libs    := $(shell root-config --libs)
glibs   := $(shell root-config --glibs)
cflags  := $(shell root-config --cflags)

#incfile := -I${PWD}/../../GemView/evio/include -Ilib/include 
#incfile := -I${PWD}/../../GemView/evio/include -I/work/prad/PRadDecoder/lib/include 
#incfile := -I${PWD}/../../GemView/evio/include -I/work/prad/xbai/PRadDecoder/lib/include 
#incfile := -I${PWD}/../../GemView/evio/include -I./PRadDecoder/lib/include 
incfile := -I/home/xbai/w/coda/common/include -I/home/xbai/w/pRad/source/PRadDecoder/lib/include -I./include

#flags   := -O3 -std=c++11 $(glibs) $(cflags) $(incfile) -L${PWD}/../../GemView/evio/lib64 -levio -levioxx -lexpat -L/work/prad/xiongw/PRadDecoder/lib -lPRadDecoder
#flags   := -O3 -std=c++11 $(glibs) $(cflags) $(incfile) -L${PWD}/../../GemView/evio/lib64 -levio -levioxx -lexpat -L/work/prad/PRadDecoder/lib -lPRadDecoder
#flags   := -O3 -std=c++11 $(glibs) $(cflags) $(incfile) -L${PWD}/../../GemView/evio/lib64 -levio -levioxx -lexpat -L/work/prad/xbai/PRadDecoder/lib -lPRadDecoder
flags   := -O3 -std=c++11 $(glibs) $(cflags) $(incfile) -L/home/xbai/w/coda/Linux-x86_64/lib -levio -levioxx -lexpat -L/home/xbai/w/pRad/source/PRadDecoder/lib -lPRadDecoder


$(Target) : $(OBJS)
	@$(cc) -o $(Target) $(OBJS) $(flags)

$(OBJDIR)/%.o: ./src/%.cxx 
	@echo Compiling $< ...
	@$(cc) -c $< -o $@ $(flags)

clean:
	@rm -f $(Target)
	@rm -f $(OBJDIR)/*.o
