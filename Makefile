################################################################################
# HiggsAnalysis/Combined Limit Makefile                                        #
#                                                                              #
# Authors: Danilo Piparo, Giovanni Petrucciani, Mingshui Chen                  #
# Revised: Nick Smith 2022                                                     #
#                                                                              #
# o Automatic compilation of new programs and classes*.                        #
# o Now generate dictionaries by genreflex                                     #
#                                                                              #
# * progs should have cpp extension, classes .cc or .cxx, and headers .h       # 
#                                                                              #
################################################################################

####  SET UP YOUR ENVIRONMENT FIRST WITH ##############################
# source env_standalone.sh
# OR
# source env_lcg.sh (if `make LCG=1` is used to build)
#######################################################################

# Hardcoded paths for standalone version identical to CMSSW 14_1_X
# These are ignored if either CONDA=1 or LCG=1 is set
BOOST = /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/boost/1.80.0-87b5de10acd2f2c8a325345ad058b814
VDT   = /cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/vdt/0.4.3-793cee1e1edef0e54b2bd5cb1f69aec9
GSL = /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/gsl/2.6-5e2ce72ea2977ff21a2344bbb52daf5c
EIGEN = /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/eigen/3bb6a48d8c171cf20b5f8e48bfb4e424fbd4f79e-3ca740c03e68b1a067f3ed0679234a78
# Compiler and flags -----------------------------------------------------------
CXX = $(shell root-config --cxx)
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs --glibs)
ROOTINC = $(shell root-config --incdir)

# CMSSW CXXFLAGS plus -Wno-unused-local-typedefs (otherwise we get a flood of messages from BOOST) plus -Wno-unused-function
CCFLAGS = -D STANDALONE $(ROOTCFLAGS) -g -fPIC -O2 -pthread -pipe -Werror=main -Werror=pointer-arith -Werror=overlength-strings -Wno-vla -Werror=overflow -ftree-vectorize -Wstrict-overflow -Werror=array-bounds -Werror=format-contains-nul -Werror=type-limits -fvisibility-inlines-hidden -fno-math-errno --param vect-max-version-for-alias-checks=50 -Xassembler --compress-debug-sections -msse3 -felide-constructors -fmessage-length=0 -Wall -Wno-non-template-friend -Wno-long-long -Wreturn-type -Wunused -Wparentheses -Wno-deprecated -Werror=return-type -Werror=missing-braces -Werror=unused-value -Werror=address -Werror=format -Werror=sign-compare -Werror=write-strings -Werror=delete-non-virtual-dtor -Werror=strict-aliasing -Werror=narrowing -Werror=unused-but-set-variable -Werror=reorder -Werror=unused-variable -Werror=conversion-null -Werror=return-local-addr -Wnon-virtual-dtor -Werror=switch -fdiagnostics-show-option -Wno-unused-local-typedefs -Wno-attributes -Wno-psabi -Wno-error=unused-variable -DBOOST_DISABLE_ASSERTS -DGNU_GCC -D_GNU_SOURCE -DBOOST_SPIRIT_THREADSAFE -DPHOENIX_THREADSAFE
LIBS = $(ROOTLIBS) -lgsl -lRooFit -lRooFitCore -lRooStats -lMinuit -lMathMore -lFoam -lHistFactory -lboost_filesystem -lboost_program_options -lboost_system -lvdt

ifeq ($(CONDA), 1)
CCFLAGS += -I${CONDA_PREFIX}/include/boost -I ${CONDA_PREFIX}/include/vdt -I ${CONDA_PREFIX}/include/gsl -I ${CONDA_PREFIX}/include/eigen3 
LIBS += -L${CONDA_PREFIX}/lib 
else ifeq ($(LCG), 1)
# for some reason, Eigen headers are nested in LCG
CCFLAGS += -I ${CPLUS_INCLUDE_PATH}/eigen3
LIBS += -L${CPLUS_INCLUDE_PATH}/../lib
else
CCFLAGS += -I$(BOOST)/include -I$(VDT)/include -I$(GSL)/include -I$(EIGEN)/include/eigen3
LIBS += -L$(BOOST)/lib -L$(VDT)/lib -L$(GSL)/lib 
endif 

# Library name -----------------------------------------------------------------
LIBNAME=HiggsAnalysisCombinedLimit
SONAME=lib$(LIBNAME).so
DICTNAME=$(LIBNAME)_xr

# Linker and flags -------------------------------------------------------------
LD = $(shell root-config --ld)
ROOTLDFLAGS   = $(shell root-config --ldflags)
ROOTLIBDIR = $(shell root-config --libdir)
# OS x specific linkage
DARWIN := $(shell uname|grep Darwin)
ifdef DARWIN
LDFLAGS       = $(ROOTLDFLAGS) -g -shared -install_name @rpath/$(SONAME) -fPIC
EXELDFLAGS    = -Wl,-rpath,'@executable_path/../lib' -Wl,-rpath,$(ROOTLIBDIR)
else
LDFLAGS       = $(ROOTLDFLAGS) -shared -Wl,-soname,$(SONAME) -Wl,-E -Wl,-z,defs -fPIC
EXELDFLAGS    = 
endif

# Directory structure ----------------------------------------------------------
PARENT_DIR = $(shell pwd)/../../
SRC_DIR = src
INC_DIR = interface
PROG_DIR = bin
SCRIPTS_DIR = scripts
PYTHON_DIR = python
# outputs
OBJ_DIR = build/obj
LIB_DIR = build/lib
EXE_DIR = build/bin


# Useful shortcuts -------------------------------------------------------------
SRCS = $(notdir $(shell ls $(SRC_DIR)/*.cc ))
SRXS = $(notdir $(shell ls $(SRC_DIR)/*.cxx ))
OBJS = $(SRCS:.cc=.o) 
OBJS += $(SRXS:.cxx=.o)
PROGS = $(notdir $(wildcard ${PROG_DIR}/*.cpp)) 
EXES = $(PROGS:.cpp=)
SCRIPTS = $(notdir $(wildcard ${SCRIPTS_DIR}/*.py)) 
PYLIB_DIR = $(LIB_DIR)/python

#Makefile Rules ---------------------------------------------------------------
.PHONY: clean exe python

all: exe python

#---------------------------------------

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)/a

$(OBJ_DIR)/a/$(DICTNAME).cc: $(SRC_DIR)/classes_def.xml | $(OBJ_DIR)
	genreflex $(SRC_DIR)/classes.h -s $< -o $@ --deep --fail_on_warnings --rootmap=$(OBJ_DIR)/a/$(DICTNAME).rootmap --rootmap-lib=$(SONAME) -Isrc -I$(PARENT_DIR)
	mv $(OBJ_DIR)/a/$(DICTNAME).rootmap $(LIB_DIR)/
	mv $(OBJ_DIR)/a/$(DICTNAME)_rdict.pcm $(LIB_DIR)/

$(OBJ_DIR)/a/%.o: $(OBJ_DIR)/a/%.cc | $(OBJ_DIR)
	$(CXX) $(CCFLAGS) -I . -I $(SRC_DIR) -I $(PARENT_DIR) -c $< -o $@

#---------------------------------------

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(INC_DIR)/%.h | $(OBJ_DIR)
	$(CXX) $(CCFLAGS) -I $(INC_DIR) -I $(SRC_DIR) -I $(PARENT_DIR) -c $< -o $@
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(SRC_DIR)/%.h | $(OBJ_DIR)
	$(CXX) $(CCFLAGS) -I $(INC_DIR) -I $(SRC_DIR) -I $(PARENT_DIR) -c $< -o $@
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx $(INC_DIR)/%.h | $(OBJ_DIR)
	$(CXX) $(CCFLAGS) -I $(INC_DIR) -I $(SRC_DIR) -I $(PARENT_DIR) -c $< -o $@

# this has no header
$(OBJ_DIR)/tdrstyle.o: $(SRC_DIR)/tdrstyle.cc
	$(CXX) $(CCFLAGS) -I $(INC_DIR) -c $< -o $@

#---------------------------------------

$(LIB_DIR):
	@mkdir -p $(LIB_DIR)

${LIB_DIR}/$(SONAME): $(addprefix $(OBJ_DIR)/,$(OBJS)) $(OBJ_DIR)/a/$(DICTNAME).o | $(LIB_DIR)
	$(LD) $(LDFLAGS) $(BOOST_INC) $^ $(SOFLAGS) -o $@ $(LIBS)

#---------------------------------------

$(EXE_DIR):
	@mkdir -p $(EXE_DIR)

exe: $(addprefix $(EXE_DIR)/,$(EXES)) $(addprefix $(EXE_DIR)/,$(SCRIPTS))
	@echo $^

$(EXE_DIR)/%: $(PROG_DIR)/%.cpp $(LIB_DIR)/$(SONAME) | $(EXE_DIR)
	@echo $<
	$(CXX) $< -o $@ $(CCFLAGS) -L $(LIB_DIR) -l $(LIBNAME) -I $(INC_DIR) -I $(SRC_DIR) -I $(PARENT_DIR) $(BOOST_INC) $(LIBS) $(EXELDFLAGS)

$(EXE_DIR)/%.py: $(SCRIPTS_DIR)/%.py | $(EXE_DIR)
	cp $< $@
# macOS System Integrity Protection unsets LD_LIBRARY_PATH for child process started by system programs
# breaking the use of /usr/bin/env in the scripts, so we hardcode the path to python executable instead
ifdef DARWIN
	sed -i "" "1s@/.*@$(shell which python)@" $@
endif

#---------------------------------------

.FORCE:

python: .FORCE | $(LIB_DIR)
	@mkdir -p $(PYLIB_DIR)/HiggsAnalysis/CombinedLimit
	@touch $(PYLIB_DIR)/__init__.py
	@touch $(PYLIB_DIR)/HiggsAnalysis/__init__.py
	@touch $(PYLIB_DIR)/HiggsAnalysis/CombinedLimit/__init__.py
	cp -r $(PYTHON_DIR)/* $(PYLIB_DIR)/HiggsAnalysis/CombinedLimit
	python3 -m compileall -q $(PYLIB_DIR)

#---------------------------------------

clean:
	@rm -rf $(OBJ_DIR) 
	@rm -rf $(EXE_DIR)
	@rm -rf $(LIB_DIR)
