# Simple Makefile to build the exampleCreateNLL helper.

CXX           ?= $(shell root-config --cxx)
CXXFLAGS      ?= -O2 -std=c++17
ROOTCFLAGS    := $(shell root-config --cflags)
ROOTLIBS      := $(shell root-config --libs)
ROOFITLIBS    := -lRooFit -lRooFitCore -lRooStats -lHistFactory

# Adjust these if your checkout or build products live elsewhere.
COMBINE_SRCDIR ?= ../..
COMBINE_LIBDIR ?= ../../build

INCLUDES := -I$(COMBINE_SRCDIR)
LIBS     := -L$(COMBINE_LIBDIR) -lHiggsAnalysisCombinedLimit
RPATH    := -Wl,-rpath,$(abspath $(COMBINE_LIBDIR))

TARGET   := exampleCreateNLL
SOURCES  := $(TARGET).cxx

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(INCLUDES) $< -o $@ $(ROOTLIBS) $(ROOFITLIBS) $(LIBS) $(RPATH)

.PHONY: clean
clean:
	rm -f $(TARGET)
