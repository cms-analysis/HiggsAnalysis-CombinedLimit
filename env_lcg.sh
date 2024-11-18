LCG_RELEASE=dev4 # or dev3 for ROOT master (dev4 is the latest ROOT stable branch)
LCG_PATH=/cvmfs/sft.cern.ch/lcg/views/$LCG_RELEASE/latest/x86_64-el9-gcc13-opt

source $LCG_PATH/setup.sh
source $LCG_PATH/bin/thisroot.sh

export PATH=$PWD/install/bin:$PATH
export LD_LIBRARY_PATH=$PWD/install/lib:$LD_LIBRARY_PATH
