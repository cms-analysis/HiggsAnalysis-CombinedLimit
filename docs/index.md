# Introduction

These pages document the
[RooStats](https://twiki.cern.ch/twiki/bin/view/RooStats/WebHome) /
[RooFit](https://root.cern.ch/roofit) - based software tools used for
statistical analysis within the [Higgs PAG](HiggsWG) - **combine**.

Combine provides a command line interface to many different statistical techniques available inside RooFit/RooStats used widely inside CMS.

The package exists in GIT under [https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit)

For more information about GIT and its usage in CMS, see [http://cms-sw.github.io/cmssw/faq.html](http://cms-sw.github.io/cmssw/faq.html)

The code can be checked out from GIT and compiled on top of a CMSSW release that includes a recent RooFit/RooStats

# Installation instructions

## Within CMSSW (recommended for CMS users)

The instructions below are for installation within a CMSSW environment. For end
users that don't need to commit or do any development, the following recipes
should be sufficient. To choose a release version, you can find the latest
releases on github under
[https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases)

### Combine v9 - recommended version

The nominal installation method is inside CMSSW. The current release targets
CMSSW `11_3_X` series because this release has both python2 and python3 ROOT
bindings, allowing a more gradual migration of user code to python3. Combine is
fully python3-compatible and can work also in 12_X releases.

```sh
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
```
Update to a recommended tag - currently the recommended tag is **v9.1.0**: [see release notes](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases/tag/v9.1.0)

```sh
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v9.1.0
scramv1 b clean; scramv1 b # always make a clean build
```

### Combine v8: `CMSSW_10_2_X` release series

Setting up the environment (once):

```sh
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
```
Update to a recommended tag - currently the recommended tag is **v8.2.0**: [see release notes](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases/tag/v8.2.0)

```sh
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.2.0
scramv1 b clean; scramv1 b # always make a clean build
```

### SLC6/CC7 release `CMSSW_8_1_X`

Setting up the environment (once):

```sh
# For CC7:
export SCRAM_ARCH=slc7_amd64_gcc530
# For SLC6:
export SCRAM_ARCH=slc6_amd64_gcc530

cmsrel CMSSW_8_1_0
cd CMSSW_8_1_0/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
```
Update to a recommended tag - currently the recommended tag for CMSSW_8_1_X is **v7.0.13**:

```sh
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v7.0.13
scramv1 b clean; scramv1 b # always make a clean build
```

## Standalone compilation

The standalone version can be easily compiled using
[cvmfs](https://cernvm.cern.ch/fs/) as it relies on dependencies which are
already installed at `/cvmfs/cms.cern.ch/`. Access to `/cvmfs/cms.cern.ch/` can
be obtained from lxplus machines or via `CernVM`. See [CernVM](CernVM.md) for
further details on the latter. In case you do not want to use the `cvmfs`
area, you will need to adapt the location of the dependencies listed in both
the `Makefile` and `env_standalone.sh` files.

```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit/ 
# git checkout <some release>
. env_standalone.sh
make -j 4
```

You will need to source `env_standalone.sh` each time you want to use the package, or add it to your login.

### Standalone compilation with LCG
LCG102 contains ROOT 6.26, which at the time of writing is not available in any CMSSW release.
This repo can be compiled against LCG102 with:
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
source env_lcg.sh 
make LCG=1 -j 8
```

### Standalone compilation with `conda`
This recipe will work both for linux and MacOS
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit

conda install --name base mamba # faster conda
mamba env create -f conda_env.yml

conda activate combine
source set_conda_env_vars.sh
# Need to reactivate
conda deactivate
conda activate combine

make CONDA=1 -j 8
```

Using combine from then on should only require sourcing the conda environment 
```
conda activate combine
```

**Note:** on OS X, `combine` can only accept workspaces, so run `text2workspace.py` first.
This is due to some ridiculous issue with child processes and `LD_LIBRARY_PATH` (see note in Makefile)

# What has changed between tags? 

You can generate a diff of any two tags (eg for `v7.0.8` and `v7.0.6`) by using following the url:

[https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/compare/v7.0.6...v7.0.7](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/compare/v7.0.6...v7.0.7)

Replace the tag names in the url to any tags you which to compare.

# For developers

We use the _Fork and Pull_ model for development: each user creates a copy of the repository on github, commits their requests there and then sends pull requests for the administrators to merge.

_Prerequisites_

1. Register on github, as needed anyway for CMSSW development: [http://cms-sw.github.io/cmssw/faq.html](http://cms-sw.github.io/cmssw/faq.html)

2. Register your SSH key on github: [https://help.github.com/articles/generating-ssh-keys](https://help.github.com/articles/generating-ssh-keys) 1 Fork the repository to create your copy of it: [https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/fork](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/fork) (more documentation at [https://help.github.com/articles/fork-a-repo](https://help.github.com/articles/fork-a-repo) )

You will now be able to browse your fork of the repository from [https://github.com/your-github-user-name/HiggsAnalysis-CombinedLimit](https://github.com/your-github-user-name/HiggsAnalysis-CombinedLimit)

# Combine Tool

An additional tool for submitting combine jobs to batch/crab, developed originally for HiggsToTauTau. Since the repository contains a certain amount of analysis-specific code, the following scripts can be used to clone it with a sparse checkout for just the core [`CombineHarvester/CombineTools`](https://github.com/cms-analysis/CombineHarvester/blob/master/CombineTools/) subpackage, speeding up the checkout and compile times:

git clone via ssh:

```sh
bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/main/CombineTools/scripts/sparse-checkout-ssh.sh)
```

git clone via https:

```sh
bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/main/CombineTools/scripts/sparse-checkout-https.sh)
```

make sure to run `scram`  to compile the `CombineTools` package.

See the [`CombineHarvester`](http://cms-analysis.github.io/CombineHarvester/) documentation pages for more details on using this tool and additional features available in the full package.


