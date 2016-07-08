Target  := main

Source  := main.cxx \
	   GEMConfigure.cxx \
	   GEMHistoManager.cxx \
	   GEMPhysHandler.cxx

OBJDIR 	:= ./obj_lib

Objs	:= $(patsubst %.cxx, %.o, $(Source))
OBJS	:= $(patsubst %.o, $(OBJDIR)/%.o, $(Objs))

cc      := g++

libs    := $(shell root-config --libs)
glibs   := $(shell root-config --glibs)
cflags  := $(shell root-config --cflags)

incfile := -I/home/xbai/w/coda/common/include -I/home/xbai/w/pRad/source/PRadEventViewer/include -I./include

flags   := -O3 -std=c++11 $(glibs) $(cflags) $(incfile) -L/home/xbai/w/coda/Linux-x86_64/lib -levio -levioxx -lexpat -L/home/xbai/w/pRad/source/PRadEventViewer/decoder/lib -lPRadDecoder


$(Target) : $(OBJS)
	@$(cc) -o $(Target) $(OBJS) $(flags)

$(OBJDIR)/%.o: ./src/%.cxx 
	@echo Compiling $< ...
	@$(cc) -c $< -o $@ $(flags)

clean:
	@rm -f $(Target)
	@rm -f $(OBJDIR)/*.o
