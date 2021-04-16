#HEADERS
CFLAGS = $(shell root-config --cflags)
LDFLAGS = $(shell root-config --ldflags)
ROOTGLIBS = $(shell root-config --glibs)
ROOTINCDIR = $(shell root-config --incdir)

EDEPSIM = /home/federico/Documents/Dune/edep-sim/edep-gcc-7-x86_64-linux-gnu

EDEPGLIBS = -L$(EDEPSIM)/lib/ -ledepsim -ledepsim_io
EDEPINCDIR = $(EDEPSIM)/include/EDepSim

#all: Digitize Reconstruct Analyze
all: Select GetDistribution
#Time MinDigitize


Select: 
	g++ src/l2g_select.cpp -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom -lEG \
	$(EDEPGLIBS) -Llib -g

GetDistribution: 
	g++ src/l2g_Distributions.cpp -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom -lEG \
	$(EDEPGLIBS) -Llib -g
	


