LCG_RELEASE=dev3 # dev3 for ROOT master (dev4 is the latest ROOT stable branch)
LCG_PATH=/cvmfs/sft.cern.ch/lcg/views/$LCG_RELEASE/latest/x86_64-el9-gcc13-opt

source $LCG_PATH/setup.sh
source $LCG_PATH/bin/thisroot.sh

export PATH=$PWD/build/bin:$PATH
export LD_LIBRARY_PATH=$PWD/build/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/build/python:$PYTHONPATH

# Run the following combine "hello world" example to test
# text2workspace.py data/tutorials/CAT23001/datacard-3-parametric-analysis.txt
# combine -M MultiDimFit -m 125.38 data/tutorials/CAT23001/datacard-3-parametric-analysis.root
