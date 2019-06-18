################################################################################
# HiggsAnalysis/Combined Limit Makefile                                        #
#                                                                              #
# Authors: Danilo Piparo, Giovanni Petrucciani, Mingshui Chen                  #
#                                                                              #
# o Automatic compilation of new programs and classes*.                        #
# o Now generate dictionaries by genreflex                                     #
#                                                                              #
# * progs should have cpp extension, classes .cc or .cxx, and headers .h       # 
#                                                                              #
################################################################################

####  SET UP YOUR ENVIRONMENT FIRST WITH ##############################
# . /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/etc/profile.d/init.sh 
# . /cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/6.02.00-odfocd5/bin/thisroot.sh
# . /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/vdt/v0.3.2-cms/etc/profile.d/init.sh 
# . /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/boost/1.51.0-cms/etc/profile.d/init.sh 
# . /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/xz/5.2.1/etc/profile.d/init.sh 
# export PATH=${PATH}:${PWD}/exe:${PWD}/scripts
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/lib
# export PYTHONPATH=${PYTHONPATH}:${PWD}/lib/python
#######################################################################

# Boost
BOOST = /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0-gnimlf
VDT   = /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/vdt/0.4.0-gnimlf
# PCRE = /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/pcre/8.37-omkpbe2
GSL = /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/gsl/2.2.1-omkpbe2
# LIBXML = /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/libxml2/2.9.1-omkpbe2/include/libxml2
# XZ = /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/xz/5.2.2-omkpbe2
# ZLIB = /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/zlib-x86_64/1.2.11-omkpbe2
# Compiler and flags -----------------------------------------------------------
CC = c++
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs --glibs)
ROOTINC = $(shell root-config --incdir)

# CCFLAGS = -D STANDALONE $(ROOTCFLAGS) -I$(BOOST)/include -I$(VDT)/include -I$(PCRE)/include -I$(GSL)/include -I$(LIBXML)/include/libxml2 -I$(XZ)/include -I$(ZLIB)/include -g -fPIC
CCFLAGS = -D STANDALONE $(ROOTCFLAGS) -I$(BOOST)/include -I$(VDT)/include -I$(GSL)/include -g -fPIC
# CMSSW CXXFLAGS plus -Wno-unused-local-typedefs (otherwise we get a flood of messages from BOOST) plus -Wno-unused-function
CCFLAGS += -O2 -pthread -pipe -Werror=main -Werror=pointer-arith -Werror=overlength-strings -Wno-vla -Werror=overflow -std=c++1z -ftree-vectorize -Wstrict-overflow -Werror=array-bounds -Werror=format-contains-nul -Werror=type-limits -fvisibility-inlines-hidden -fno-math-errno --param vect-max-version-for-alias-checks=50 -Xassembler --compress-debug-sections -msse3 -felide-constructors -fmessage-length=0 -Wall -Wno-non-template-friend -Wno-long-long -Wreturn-type -Wunused -Wparentheses -Wno-deprecated -Werror=return-type -Werror=missing-braces -Werror=unused-value -Werror=address -Werror=format -Werror=sign-compare -Werror=write-strings -Werror=delete-non-virtual-dtor -Werror=strict-aliasing -Werror=narrowing -Werror=unused-but-set-variable -Werror=reorder -Werror=unused-variable -Werror=conversion-null -Werror=return-local-addr -Wnon-virtual-dtor -Werror=switch -fdiagnostics-show-option -Wno-unused-local-typedefs -Wno-attributes -Wno-psabi -Wno-error=unused-variable -DBOOST_DISABLE_ASSERTS -DGNU_GCC -D_GNU_SOURCE -DBOOST_SPIRIT_THREADSAFE -DPHOENIX_THREADSAFE
LIBS = $(ROOTLIBS) -L$(BOOST)/lib -L$(VDT)/lib -L$(GSL)/lib -lgsl -l RooFit -lRooFitCore -l RooStats -l Minuit -lMathMore -l Foam -lHistFactory -lboost_filesystem -lboost_program_options -lboost_system -lvdt

# Library name -----------------------------------------------------------------
LIBNAME=HiggsAnalysisCombinedLimit
SONAME=lib$(LIBNAME).so
DICTNAME=$(LIBNAME)_xr

# Linker and flags -------------------------------------------------------------
LD = g++
ROOTLDFLAGS   = $(shell root-config --ldflags)
LDFLAGS       = $(ROOTLDFLAGS) -shared -Wl,-soname,$(SONAME) -Wl,-E -Wl,-z,defs -fPIC

# Directory structure ----------------------------------------------------------
PARENT_DIR = $(shell pwd)/../../
SRC_DIR = src
INC_DIR = interface
LIB_DIR = lib
PROG_DIR = bin
EXE_DIR = exe
OBJ_DIR = obj


# Useful shortcuts -------------------------------------------------------------
SRCS = $(notdir $(shell ls $(SRC_DIR)/*.cc ))
SRXS = $(notdir $(shell ls $(SRC_DIR)/*.cxx ))
OBJS = $(SRCS:.cc=.o) 
OBJS += $(SRXS:.cxx=.o)
PROGS = $(notdir $(wildcard ${PROG_DIR}/*.cpp)) 
EXES = $(PROGS:.cpp=)

#Makefile Rules ---------------------------------------------------------------
.PHONY: clean dirs dict obj lib exe debug


all: dirs dict obj lib exe compile_python

#---------------------------------------

dirs:
	@mkdir -p $(OBJ_DIR)/a
	@mkdir -p $(SRC_DIR)
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(EXE_DIR)
	@mkdir -p $(LIB_DIR)/python/HiggsAnalysis
	@ln -sd ../../../python $(LIB_DIR)/python/HiggsAnalysis/CombinedLimit || /bin/true
	@touch $(LIB_DIR)/python/__init__.py
	@touch $(LIB_DIR)/python/HiggsAnalysis/__init__.py
	@touch $(LIB_DIR)/python/HiggsAnalysis/CombinedLimit/__init__.py

#---------------------------------------

dict: dirs $(OBJ_DIR)/a/$(DICTNAME).cc
$(OBJ_DIR)/a/$(DICTNAME).cc : $(SRC_DIR)/classes_def.xml
# 	@echo "\n*** Generating dictionaries ..."
	genreflex $(SRC_DIR)/classes.h -s $(SRC_DIR)/classes_def.xml -o $(OBJ_DIR)/a/$(DICTNAME).cc --deep --fail_on_warnings --rootmap=$(OBJ_DIR)/a/$(DICTNAME).rootmap --rootmap-lib=$(SONAME) -I$(PARENT_DIR) 
	mv $(OBJ_DIR)/a/$(DICTNAME).rootmap $(LIB_DIR)/
	mv $(OBJ_DIR)/a/$(DICTNAME)_rdict.pcm $(LIB_DIR)/

#---------------------------------------

obj: dict 
# 	@echo "\n*** Compiling ..."
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cc $(INC_DIR)/%.h
	$(CC) $(CCFLAGS) -I $(INC_DIR) -I $(SRC_DIR) -I $(PARENT_DIR) -c $< -o $@
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cc $(SRC_DIR)/%.h
	$(CC) $(CCFLAGS) -I $(INC_DIR) -I $(SRC_DIR) -I $(PARENT_DIR) -c $< -o $@
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cxx $(INC_DIR)/%.h
	$(CC) $(CCFLAGS) -I $(INC_DIR) -I $(SRC_DIR) -I $(PARENT_DIR) -c $< -o $@
$(OBJ_DIR)/a/%.o : $(OBJ_DIR)/a/%.cc 
	$(CC) $(CCFLAGS) -I . -I $(SRC_DIR) -I $(PARENT_DIR) -c $< -o $@


# this has no header
$(OBJ_DIR)/tdrstyle.o: $(SRC_DIR)/tdrstyle.cc
	$(CC) $(CCFLAGS) -I $(INC_DIR) -c $< -o $@

#---------------------------------------

lib: dirs ${LIB_DIR}/$(SONAME)
${LIB_DIR}/$(SONAME):$(addprefix $(OBJ_DIR)/,$(OBJS)) $(OBJ_DIR)/a/$(DICTNAME).o 
#	@echo "\n*** Building $(SONAME) library:"
	$(LD) $(LDFLAGS) $(BOOST_INC) $(addprefix $(OBJ_DIR)/,$(OBJS)) $(OBJ_DIR)/a/$(DICTNAME).o $(SOFLAGS) -o $@ $(LIBS)

#---------------------------------------

exe: $(addprefix $(EXE_DIR)/,$(EXES))
# 	@echo "\n*** Compiling executables ..."
$(EXE_DIR)/% : $(PROG_DIR)/%.cpp lib
	$(CC) $< -o $@ $(CCFLAGS) -L $(LIB_DIR) -l $(LIBNAME) -I $(INC_DIR) -I $(SRC_DIR) -I $(PARENT_DIR) $(BOOST_INC) $(LIBS)

compile_python:
	@python -m compileall -q python 

#---------------------------------------

clean:
# 	@echo "*** Cleaning all directories and dictionaries ..."
	@rm -rf $(OBJ_DIR) 
	@rm -rf $(EXE_DIR)
	@rm -rf $(LIB_DIR)
	@rm -rf python/*pyc python/*/*pyc

#---------------------------------------

debug:
	@echo "OBJS: $(OBJS)"
	@echo "SRCS: $(SRCS)"
