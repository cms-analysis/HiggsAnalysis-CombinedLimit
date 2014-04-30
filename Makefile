################################################################################
# HiggsAnalysis/Combined Limit Makefile                                        #
#                                                                              #
# Authors: Danilo Piparo, Giovanni Petrucciani                                 #
#                                                                              #
# o Automatic compilation of new programs and classes*.                        #
# o Automatic generation of CINT dictionaries via rootcint.                    #
#                                                                              #
# * progs should have cpp extension, classes .cc or .cxx, and headers .h       # 
#                                                                              #
################################################################################

####  SET UP YOUR ENVIRONMENT FIRST WITH ##############################
# . /afs/cern.ch/cms/slc6_amd64_gcc481/external/gcc/4.8.1/etc/profile.d/init.sh 
# . /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.17/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh 
# . /afs/cern.ch/cms/slc6_amd64_gcc481/cms/vdt/v0.3.2-cms/etc/profile.d/init.sh 
# . /afs/cern.ch/cms/slc6_amd64_gcc481/external/boost/1.51.0-cms/etc/profile.d/init.sh 
# export PATH=${PATH}:${PWD}/exe:${PWD}/scripts
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/lib
# export PYTHONPATH=${PYTHONPATH}:${PWD}/lib/python
#######################################################################

# Boost
BOOST = /afs/cern.ch/cms/slc6_amd64_gcc481/external/boost/1.51.0-cms
VDT   = /afs/cern.ch/cms/slc6_amd64_gcc481/cms/vdt/v0.3.2-cms

# Compiler and flags -----------------------------------------------------------
CC = g++

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs --glibs)
ROOTINC = $(shell root-config --incdir)

CCFLAGS = -D STANDALONE $(ROOTCFLAGS) -I$(BOOST)/include -I$(VDT)/include -g -fPIC
# CMSSW CXXFLAGS plus -Wno-unused-local-typedefs (otherwise we get a flood of messages from BOOST) plus -Wno-unused-function
CCFLAGS += -O2 -pedantic -pthread -pipe -Wno-vla -Werror=overflow -Wstrict-overflow -std=c++0x -msse3 -ftree-vectorize -Wno-strict-overflow -Werror=array-bounds -Werror=format-contains-nul -Werror=type-limits -fvisibility-inlines-hidden -fno-math-errno --param vect-max-version-for-alias-checks=50 -fipa-pta -felide-constructors -fmessage-length=0 -ftemplate-depth-300 -Wall -Wno-non-template-friend -Wno-long-long -Wreturn-type -Wunused -Wparentheses -Wno-deprecated -Werror=return-type -Werror=missing-braces -Werror=unused-value -Werror=address -Werror=format -Werror=sign-compare -Werror=write-strings -Werror=delete-non-virtual-dtor -Werror=maybe-uninitialized -Werror=strict-aliasing -Werror=narrowing -Werror=uninitialized -Werror=unused-but-set-variable -Werror=reorder -Werror=unused-variable -Werror=conversion-null -Werror=switch -fdiagnostics-show-option -DBOOST_DISABLE_ASSERTS -Wno-unused-local-typedefs -Wno-unused-function
LIBS = $(ROOTLIBS) -L$(BOOST)/lib -L$(VDT)/lib -l RooFit -lRooFitCore -l RooStats -l Minuit -l Foam -lHistFactory -lboost_filesystem -lboost_program_options -lboost_system -lvdt

# Library name -----------------------------------------------------------------
LIBNAME=HiggsAnalysisCombinedLimit
SONAME=lib$(LIBNAME).so

# Linker and flags -------------------------------------------------------------
LD = g++
ROOTLDFLAGS   = $(shell root-config --ldflags)
LDFLAGS       = $(ROOTLDFLAGS) -rdynamic -shared -Wl,-soname,$(SONAME) -fPIC

# Dictionaries filename --------------------------------------------------------
DICTNAME=cintdictionary

# Directory structure ----------------------------------------------------------
SRC_DIR = src
INC_DIR = interface
LIB_DIR = lib
PROG_DIR = bin
EXE_DIR = exe
OBJ_DIR = obj


# Useful shortcuts -------------------------------------------------------------
SRCS = $(notdir $(shell ls $(SRC_DIR)/*.cc|grep -v $(DICTNAME) ))
SRXS = $(notdir $(shell ls $(SRC_DIR)/*.cxx|grep -v $(DICTNAME) ))
SRCS += $(DICTNAME).cc
OBJS = $(SRCS:.cc=.o) 
OBJS += $(SRXS:.cxx=.o)
PROGS = $(notdir $(wildcard ${PROG_DIR}/*.cpp)) 
EXES = $(PROGS:.cpp=)

# Classes with dicts -----------------------------------------------------------
DICTHDRS = $(notdir $(shell grep -l ClassDef interface/*h ))
# Classes with no ClassDef but still with a dict ------------------------------
DICTHDRS += SequentialMinimizer.h
# Functions with dicts ---------------------------------------------------------
DICTHDRS += th1fmorph.h

#Makefile Rules ---------------------------------------------------------------
.PHONY: clean dirs dict obj lib exe debug


all: dirs dict obj lib exe

#---------------------------------------

dirs:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(EXE_DIR)
	@mkdir -p $(LIB_DIR)/python/HiggsAnalysis
	@ln -sd ../../../python $(LIB_DIR)/python/HiggsAnalysis/CombinedLimit || /bin/true
	@touch $(LIB_DIR)/python/__init__.py
	@touch $(LIB_DIR)/python/HiggsAnalysis/__init__.py
	@touch $(LIB_DIR)/python/HiggsAnalysis/CombinedLimit/__init__.py

#---------------------------------------

dict: dirs $(SRC_DIR)/$(DICTNAME).cc
$(SRC_DIR)/$(DICTNAME).cc : $(SRC_DIR)/CombinedLimit_LinkDef.h 
# 	@echo "\n*** Generating dictionaries ..."
	rootcint -f $(SRC_DIR)/$(DICTNAME).cc -c -p -I$(INC_DIR) -I$(SRC_DIR) -I$(ROOTINC) $(DICTHDRS) $(SRC_DIR)/CombinedLimit_LinkDef.h
	mv $(SRC_DIR)/$(DICTNAME).h $(INC_DIR)/$(DICTNAME).h 

#---------------------------------------

obj: dict 
# 	@echo "\n*** Compiling ..."
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cc $(INC_DIR)/%.h
	$(CC) $(CCFLAGS) -I $(INC_DIR) -I $(SRC_DIR) -c $< -o $@
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cc $(SRC_DIR)/%.h
	$(CC) $(CCFLAGS) -I $(INC_DIR) -I $(SRC_DIR) -c $< -o $@
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cxx $(INC_DIR)/%.h
	$(CC) $(CCFLAGS) -I $(INC_DIR) -I $(SRC_DIR) -c $< -o $@

# this has no header
$(OBJ_DIR)/tdrstyle.o: $(SRC_DIR)/tdrstyle.cc
	$(CC) $(CCFLAGS) -I $(INC_DIR) -c $< -o $@

#---------------------------------------

lib: dirs ${LIB_DIR}/$(SONAME)
${LIB_DIR}/$(SONAME):$(addprefix $(OBJ_DIR)/,$(OBJS)) 
# 		@echo "\n*** Building $(SONAME) library:"
		$(LD) $(LDFLAGS) $(BOOST_INC) $(addprefix $(OBJ_DIR)/,$(OBJS))  $(SOFLAGS) -o $@ $(LIBS)

#---------------------------------------

exe: $(addprefix $(EXE_DIR)/,$(EXES))
# 	@echo "\n*** Compiling executables ..."
$(addprefix $(EXE_DIR)/,$(EXES)) : $(addprefix $(PROG_DIR)/,$(PROGS))
	$(CC) $< -o $@ $(CCFLAGS) -L $(LIB_DIR) -l $(LIBNAME) -I $(INC_DIR) $(BOOST_INC) $(LIBS)


#---------------------------------------

clean:
# 	@echo "*** Cleaning all directories and dictionaries ..."
	@rm -rf $(OBJ_DIR) 
	@rm -rf $(EXE_DIR)
	@rm -rf $(LIB_DIR)
	@rm -rf $(SRC_DIR)/$(DICTNAME).cc
	@rm -rf $(INC_DIR)/$(DICTNAME).h

#---------------------------------------

debug:
	@echo "OBJS: $(OBJS)"
	@echo "SRCS: $(SRCS)"
