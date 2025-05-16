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

### Combine v10 - recommended version

The nominal installation method is inside CMSSW. The current release targets
the CMSSW `14_1_X` series because of the recent switch to el9 at lxplus machines.



```sh
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
```
Update to a recommended tag - currently the recommended tag is **v10.2.0**: [see release notes](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases/tag/v10.2.0)

```sh
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v10.2.0
scramv1 b clean; scramv1 b # always make a clean build
```

### Combine v9 

The nominal installation method is inside CMSSW. The current release targets
the CMSSW `11_3_X` series because this release has both python2 and python3 ROOT
bindings, allowing a more gradual migration of user code to python3. <span style="font-variant:small-caps;">Combine</span> is
fully python3-compatible and, with some adaptations, can also work in 12_X releases. 

CMSSW `11_3_X` runs on slc7, which can be setup using apptainer ([see detailed instructions](http://cms-sw.github.io/singularity.html)):
```sh
cmssw-el7
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
```
Update to a recommended tag - currently the recommended tag is **v9.2.1**: [see release notes](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases/tag/v9.2.1)

```sh
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v9.2.1
scramv1 b clean; scramv1 b # always make a clean build
```

#### Combine v8: `CMSSW_10_2_X` release series

Setting up the environment (once):

```sh
cmssw-el7
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

#### SLC6/CC7 release `CMSSW_8_1_X`

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
scramv1 b clean; scramv1 b # always make a clean build
```

### Oustide of CMSSW (recommended for non-CMS users)

Pre-compiled versions of the tool are available as container images from the [CMS cloud](https://gitlab.cern.ch/cms-cloud/combine-standalone/container_registry/15235). These containers can be downloaded and run using [Docker](https://cms-opendata-guide.web.cern.ch/tools/docker/). If you have docker running you can pull and run the image using, 

```sh
docker run --name combine -it gitlab-registry.cern.ch/cms-cloud/combine-standalone:<tag>
```
where you must replace `<tag>` with a particular version of the tool - eg - `v9.2.1`. See the top of this page for the latest recommended versions. 

You will now have the compiled <span style="font-variant:small-caps;">Combine</span> binary available as well as the complete package of tool. 
The container can be re-started using `docker start -i combine`. 

#### Standalone compilation

The standalone version can be easily compiled using
[cvmfs](https://cernvm.cern.ch/fs/) as it relies on dependencies that are
already installed at `/cvmfs/cms.cern.ch/`. Access to `/cvmfs/cms.cern.ch/` can
be obtained from lxplus machines or via `CernVM`. See [CernVM](CernVM.md) for
further details on the latter. In case you do not want to use the `cvmfs`
area, you will need to adapt the locations of the dependencies listed in both
the `Makefile` and `env_standalone.sh` files.

```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit/ 
# git checkout <some release>
. env_standalone.sh
make -j 4
```

You will need to source `env_standalone.sh` each time you want to use the package, or add it to your login environment.

##### Compilation of slc7 compatible versions

For <span style="font-variant:small-caps;">Combine</span> versions before v10 release you will need to do the compilation in an slc7 environment using apptainer. You can then source the standalone script outside of the apptainer.
On lxplus this can be done as follows:

```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit/ 
# git checkout <some release>
cmssw-el7
. env_standalone.sh
make -j 4
exit
source . env_standalone.sh
```

##### Standalone compilation with LCG
For compilation outside of CMSSW, for example to use ROOT versions not yet available in CMSSW, one can compile against LCG releases. The current default is to compile with LCG_102, which contains ROOT 6.26:
```sh
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
source env_lcg.sh 
make LCG=1 -j 8
```
To change the LCG version, edit `env_lcg.sh`. 

The resulting binaries can be moved for use in a
batch job if the following files are included in the job tarball:
```sh
tar -zcf Combine_LCG_env.tar.gz build interface src/classes.h --exclude=obj
```

##### Standalone compilation with `conda`
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

Using <span style="font-variant:small-caps;">Combine</span> from then on should only require sourcing the conda environment 
```
conda activate combine
```

**Note:** on OS X, <span style="font-variant:small-caps;">Combine</span> can only accept workspaces, so run `text2workspace.py` first.
This is due to an issue with child processes and `LD_LIBRARY_PATH` (see note in Makefile)

##### Standalone compilation with CernVM 

<span style="font-variant:small-caps;">Combine</span>, either standalone or not, can be compiled via CVMFS using access to `/cvmfs/cms.cern.ch/`  obtained using a virtual machine - [`CernVM`](https://cernvm.cern.ch/). To use `CernVM` You should have access to CERN IT resources. If you are a CERN user you can use your account, otherwise you can request a lightweight account.
If you have a CERN user account, we strongly suggest you simply run one of the other standalone installations, which are simpler and faster than using a VM.

You should have a working VM on your local machine, compatible with CernVM, such as `VirtualBox`. All the required software can be downloaded [here](https://cernvm.cern.ch/appliance/).
At least 2GB of disk space should be reserved on the virtual machine for <span style="font-variant:small-caps;">Combine</span> to work properly and the machine must be contextualized to add the `CMS` group to CVMFS. A minimal working setup is described below.

0. Download the CernVM-launcher for your operating system, following the instructions available [`here`] for your operating system (https://cernvm.readthedocs.io/en/stable/cpt-launch.html#installation

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
scram b
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
