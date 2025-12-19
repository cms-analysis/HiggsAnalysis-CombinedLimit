hide:
    - navigation

# Introduction

<img align="left" src="logo.png">
These pages document the
[RooStats](https://twiki.cern.ch/twiki/bin/view/RooStats/WebHome) /
[RooFit](https://root.cern.ch/roofit) - based software tool used for
statistical analysis within the CMS experiment - <span style="font-variant:small-caps;">Combine</span>. Note that while this tool was originally developed in the Higgs Physics Analysis Group (PAG), its usage is now widespread within CMS.

<span style="font-variant:small-caps;">Combine</span> provides a command-line interface to many different statistical techniques, available inside RooFit/RooStats, that are used widely inside CMS.

The package exists on GitHub under [https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit)

For more information about Git, GitHub and its usage in CMS, see [http://cms-sw.github.io/cmssw/faq.html](http://cms-sw.github.io/cmssw/faq.html)

The code can be checked out from GitHub and compiled on top of a CMSSW release that includes a recent RooFit/RooStats, or via standalone compilation without CMSSW dependencies. See the instructions for installation of <span style="font-variant:small-caps;">Combine</span> below.


## Installation instructions

Installation instructions and recommended versions can be found below.
Since v9.0.0, the versioning follows the [semantic versioning 2.0.0 standard](https://semver.org/).
Earlier versions are not guaranteed to follow the standard.

### Within CMSSW (recommended for CMS users)

The instructions below are for installation within a CMSSW environment. For end
users that do not need to commit or do any development, the following recipes
should be sufficient. To choose a release version, you can find the latest
releases on github under
[https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases)

#### Combine v10 - recommended version

The nominal installation method is inside CMSSW. The current release targets
the CMSSW `14_1_X` series because of the recent switch to el9 at lxplus machines.

Currently, the recommended tag is **v10.4.2**: [see release notes](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases/tag/v10.4.2)
The `git clone` command below contains this tag and is optimised to reduce disk usage.

```sh
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git -c advice.detachedHead=false clone --depth 1 --branch v10.4.2 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
scramv1 b clean; scramv1 b -j$(nproc --ignore=2) # always make a clean build, with n - 2 cores on the system
```

#### Legacy versions

##### Combine v9

The nominal installation method is inside CMSSW. The current release targets
the CMSSW `11_3_X` series because this release has both python2 and python3 ROOT
bindings, allowing a more gradual migration of user code to python3. <span style="font-variant:small-caps;">Combine</span> is
fully python3-compatible and, with some adaptations, can also work in 12_X releases.

CMSSW `11_3_X` runs on slc7, which can be setup using apptainer ([see detailed instructions](http://cms-sw.github.io/singularity.html)).
Currently, the recommended tag is **v9.2.1**: [see release notes](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases/tag/v9.2.1)
The `git clone` command below contains this tag and is optimised to reduce disk usage.

```sh
cmssw-el7
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
git -c advice.detachedHead=false clone --depth 1 --branch v9.2.1 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
scramv1 b clean; scramv1 b -j$(nproc --ignore=2) # always make a clean build, with n - 2 cores on the system
```

##### Combine v8: `CMSSW_10_2_X` release series

Setting up the environment (once) is described below.
Currently, the recommended tag is **v8.2.0**: [see release notes](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases/tag/v8.2.0)
The `git clone` command below contains this tag and is optimised to reduce disk usage.

```sh
cmssw-el7
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git -c advice.detachedHead=false clone --depth 1 --branch v8.2.0 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
scramv1 b clean; scramv1 b -j$(nproc --ignore=2) # always make a clean build, with n - 2 cores on the system
```

##### SLC6/CC7 release `CMSSW_8_1_X`

Setting up OS using apptainer ([see detailed instructions](http://cms-sw.github.io/singularity.html)):

```sh
# For CC7:
cmssw-el7
# For SLC6:
cmssw-el6

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
scramv1 b clean; scramv1 b -j$(nproc --ignore=2) # always make a clean build, with n - 2 cores on the system
```

### Outside of CMSSW

#### Standalone compilation

Combine can be built as a CMake project. It has four required dependencies on the C++ side:

  1. [ROOT](https://root.cern), mostly for RooFit
  2. The [boost](https://www.boost.org) library, for command line options parsing
  3. [Eigen](https://eigen.tuxfamily.org) for linear algebra, used in the `CMSInterferenceFunc` and `RooSplineND`

There are two more optional dependencies:

  1. The [vdt](https://github.com/dpiparo/vdt) library for fast vectorized math (can be disabled with the CMake configuration option `-DUSE_VDT=FALSE`)
  2. The [gtest](https://github.com/google/googletest) library for unit tests, if you build with `-DBUILD_TESTS=TRUE`

Any environment that provides the dependencies can be used to build combine.

To build, run the following commands inside the cloned repository:
```bash
mkdir build
cd build
cmake .. # additional CMake configuration options like -DUSE_VDT=FALSE go here
cmake --build . -j8
cd ..
```

To use your build of Combine, you have to append to the following environment variables:
```bash
export PATH=$PWD/build/bin:$PATH
export LD_LIBRARY_PATH=$PWD/build/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/build/python:$PYTHONPATH
```

For advanced users or packagers who want to install the build, there are some more relevant options to steer the CMake installation step:

* `CMAKE_INSTALL_BINDIR`: the binary directory inside the install prefix for `combine` and Python scripts like `text2workspace.py`  (`bin/` by default)
* `CMAKE_INSTALL_LIBDIR`: shared library directory (`lib/` by default)
* `CMAKE_INSTALL_PYTHONDIR`: python module directory (`python/` by default)
* `CMAKE_INSTALL_INCLUDEDIR`: header file directory (`include/` by default)

##### Standalone compilation with LCG

A typical environment that can be used on `lxplus` is the [LCG software stack](https://lcginfo.cern.ch/).

It can be activated as follows
```bash
LCG_RELEASE=LCG_106 # includes ROOT 6.32, like CMSSW_14_1_0_pre4
# LCG_RELEASE=dev3/latest # includes nightly build of ROOT master, useful for development
LCG_PATH=/cvmfs/sft.cern.ch/lcg/views/$LCG_RELEASE/x86_64-el9-gcc13-opt

source $LCG_PATH/setup.sh
source $LCG_PATH/bin/thisroot.sh
```

After activating the environment, you can follow the usual CMake build procedure explained above.

##### Standalone compilation with `conda` (CMake-based)
This recipe mirrors the setup used in our GitHub Actions builds and works on both Linux and macOS (Intel or Apple silicon):

```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit

# configure conda-forge as the preferred channel
conda config --set channel_priority strict
conda config --add channels conda-forge

# create and activate the environment
conda create -n combine python=3.12 root=6.34 gsl boost-cpp vdt eigen tbb cmake ninja
conda activate combine

# configure and build with CMake
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_INSTALL_PYTHONDIR=lib/python3.12/site-packages -DUSE_VDT=OFF
cmake --build build -j$(nproc --ignore=2)
cmake --install build
```

After installation the binaries and Python modules live inside the environment, so a new shell only requires:

```
conda activate combine
```

#### Pre-compiled image with CVMFS dependency 

Pre-compiled versions of the tool are available as container images from the [GitLab CMS-analysis repository](http://gitlab-registry.cern.ch/cms-analysis/general/combine-container). This container is built together with the [`CombineHarvester`](http://cms-analysis.github.io/CombineHarvester/) package and is **recommended for users with CVMFS access**. 

```shell
export APPTAINER_CACHEDIR="/tmp/$(whoami)/apptainer_cache"
apptainer shell -B /cvmfs -B /eos -B /afs /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/combine-container:latest /bin/bash
```

Then in the Apptainer shell:

```shell
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home/cmsusr/CMSSW_14_1_0_pre4/
cmsenv  # Ignore errors
```

Alternatively, you can set up Combine with a single-line command which relies on the same image:

```
source /cvmfs/cms.cern.ch/cat/combine_env.sh
```
Note that the batch submission is not supported with this setup; if needed refer to the [default setup instructions](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#within-cmssw-recommended-for-cms-users). 

#### Pre-compiled versions in Docker

Pre-compiled versions of the tool are available as container images from the [CMS cloud](https://gitlab.cern.ch/cms-cloud/combine-standalone/container_registry/15235). These containers can be downloaded and run using [Docker](https://cms-opendata-guide.web.cern.ch/tools/docker/). If you have docker running you can pull and run the image using,

```sh
docker run --name combine -it gitlab-registry.cern.ch/cms-cloud/combine-standalone:<tag>
```
where you must replace `<tag>` with a particular version of the tool. Containers are available from `v9.2.1` and `v10.0.1` onwards. If no tag is specified the latest version of the container will be loaded, which is `v10.4.2` at the moment.

You will now have the compiled <span style="font-variant:small-caps;">Combine</span> binary available as well as the complete package of tool.
The container can be re-started using `docker start -i combine`.

#### Standalone compilation with CernVM

<span style="font-variant:small-caps;">Combine</span>, either standalone or not, can be compiled via CVMFS using access to `/cvmfs/cms.cern.ch/`  obtained using a virtual machine - [`CernVM`](https://cernvm.cern.ch/). To use `CernVM` You should have access to CERN IT resources. If you are a CERN user you can use your account, otherwise you can request a lightweight account.
If you have a CERN user account, we strongly suggest you simply run one of the other standalone installations, which are simpler and faster than using a VM.

You should have a working VM on your local machine, compatible with CernVM, such as `VirtualBox`. All the required software can be downloaded [here](https://cernvm.cern.ch/appliance/).
At least 2GB of disk space should be reserved on the virtual machine for <span style="font-variant:small-caps;">Combine</span> to work properly and the machine must be contextualized to add the `CMS` group to CVMFS. A minimal working setup is described below.

0. Download the CernVM-launcher for your operating system, following the instructions available [here](https://cernvm.readthedocs.io/en/stable/cpt-launch.html#installation) for your operating system

1. Prepare a CMS context. You can use the CMS open data one already available on gitHub:
```wget https://raw.githubusercontent.com/cernvm/public-contexts/master/cms-opendata-2011.context)```

2. Launch the virtual machine ```cernvm-launch create --name combine --cpus 2 cms-opendata-2011.context```

3. In the VM, proceed with an installation of combine

Installation through CernVM is maintained on a best-effort basis and these instructions may not be up to date.

## What has changed between tags?

You can generate a diff of any two tags (eg for `v9.2.1` and `v9.2.0`) by using the following url:

[https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/compare/v9.2.0...v9.2.1](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/compare/v9.2.0...v9.2.1)

Replace the tag names in the url to any tags you would like to compare.

## For developers

We use the _Fork and Pull_ model for development: each user creates a copy of the repository on GitHub, commits their requests there, and then sends pull requests for the administrators to merge.

_Prerequisites_

1. Register on GitHub, as needed anyway for CMSSW development: [http://cms-sw.github.io/cmssw/faq.html](http://cms-sw.github.io/cmssw/faq.html)

2. Register your SSH key on GitHub: [https://help.github.com/articles/generating-ssh-keys](https://help.github.com/articles/generating-ssh-keys)

3. Fork the repository to create your copy of it: [https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/fork](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/fork) (more documentation at [https://help.github.com/articles/fork-a-repo](https://help.github.com/articles/fork-a-repo) )

You will now be able to browse your fork of the repository from [https://github.com/your-github-user-name/HiggsAnalysis-CombinedLimit](https://github.com/your-github-user-name/HiggsAnalysis-CombinedLimit)

We strongly encourage you to contribute any developments you make back to the main repository.
See [contributing.md](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/contributing.md) for details about contributing.

## CombineHarvester/CombineTools

CombineHarvester/CombineTools is a package for the creation of datacards/workspaces used with <span style="font-variant:small-caps;">Combine v10</span> for a number of analyses in CMS. See the [`CombineHarvester`](http://cms-analysis.github.io/CombineHarvester/) documentation pages for more details on using this tool and additional features available in the full package.

This package also comes with useful features for <span style="font-variant:small-caps;">Combine</span> such as the automated datacard validation (see [instructions](docs/part3/validation)). The repository can be checked out and compiled using,

```sh
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
scram b -j$(nproc --ignore=2)
```

See the [`CombineHarvester`](http://cms-analysis.github.io/CombineHarvester/) documentation for full instructions and reccomended versions.

!!! info
    Starting with <span style="font-variant:small-caps;">Combine v10</span>, specific ombineTool functionalities for job submition and parallelization (`combineTool.py`) as well as many plotting functions have been integrated into the <span style="font-variant:small-caps;">Combine</span> package. For these tasks you no longer have to follow the instructions above.


## Citation

If you use <span style="font-variant:small-caps;">Combine</span>, please cite the following CMS publication [here](https://arxiv.org/abs/2404.06614).

<details>
<summary><b>Show BibTex Entry</b></summary>
```
@article{
    CMS:2024onh,
    author = "Hayrapetyan, Aram and others",
    collaboration = "CMS",
    title = "The {CMS} statistical analysis and combination tool: {\textsc{Combine}}",
    eprint = "2404.06614",
    archivePrefix = "arXiv",
    primaryClass = "physics.data-an",
    reportNumber = "CMS-CAT-23-001, CERN-EP-2024-078",
    year = "2024",
    journal = "Comput. Softw. Big Sci.",
    doi = "10.1007/s41781-024-00121-4",
    volume = "8",
    pages = "19"
}
```
</details>
