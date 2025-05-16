#!/usr/bin/env bash

# Need to run with "./demo_lcg.sh" on lxplus9

# Setup LGC release and paths
source ./env_lcg.sh

# For the python part, it is recommended to install it as a python package
# into your virtual environment (see setup.py). Do this AFTER the LCG
# environment setup to use a consistent Python version.

mkdir -p .python/HiggsAnalysis
ln -s $PWD/python $PWD/.python/HiggsAnalysis/CombinedLimit
touch $PWD/.python/HiggsAnalysis/__init__.py
touch $PWD/.python/HiggsAnalysis/CombinedLimit/__init__.py

python3 -m venv virtualenv
source virtualenv/bin/activate

python3 -m pip install -e .

# Now, build the C++ part of combine
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DINSTALL_PYTHON=FALSE ..
make install -j4
cd ..

# Run the following combine "hello world" example
text2workspace.py data/tutorials/CAT23001/datacard-3-parametric-analysis.txt
combine -M MultiDimFit -m 125.38 data/tutorials/CAT23001/datacard-3-parametric-analysis.root
