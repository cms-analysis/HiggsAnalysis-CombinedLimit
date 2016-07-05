Combine Documentation {#mainpage}
=================================

Introduction
------------

The package now exists in GIT under <https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit>
Previously it was hosted on CVS under [HiggsAnalysis/CombinedLimit](http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/HiggsAnalysis/CombinedLimit/).

**For more information about GIT and its usage in CMS, see http://cms-sw.github.io/cmssw/faq.html**

The code can be checked out from GIT and compiled on top of a CMSSW release that includes a recent RooFit / [RooStats](RooStats.WebHome) (see below).

Setting up the environment
--------------------------

This package requires Scientific Linux 5 (SL5); as you cannot log anymore into lxplus5.cern.ch, you will need to create a virtual machine (VM) with SL5 installed to perform the operations below. You can find instructions on how to create a VM using the link [here](CMSAISL5InstallationAndConfiguration). With the VM you will be able to log into your afs area, and set up the environemnt. If git does not work under the VM, then you can log into lxplus, get the package through git, then go back on the virtual machine to compile it.

### GIT recipe (the only supported recipe now)

#### %BLUE%For end-users that don't need to commit or do any development<span class="twiki-macro ENDCOLOR"></span>

This is the simplest recipe and does not need a github account, nor to know anything about git besides these commands below.

You can find the latest releases on github under <https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases>

##### ROOT6 SLC6 release (CMSSW\_7\_4\_X)

**Setting up the environment (once)**

    export SCRAM_ARCH=slc6_amd64_gcc491
    cmsrel CMSSW_7_4_7
    cd CMSSW_7_4_7/src 
    cmsenv
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit

**Update to a reccomended tag - currently the reccomended tag is v6.2.1**

    cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
    git fetch origin
    git checkout v6.2.1
    scramv1 b clean; scramv1 b # always make a clean build, as scram doesn't always see updates to src/LinkDef.h

-   when running with the HybridNew calculator, you can safely ignore the following warning from RooFit
    **`WARNING:Eval -- RooStatsUtils::MakeNuisancePdf - no constraints found on nuisance parameters in the input model`**

**Combine Tool**

-   An additional tool for submitting combine jobs to batch/crab, developed originally for HiggsToTauTau. Since the repository contains a certain amount of analysis-specific code, the following scripts can be used to clone it with a sparse checkout for just the core CombineHarvester/CombineTools subpackage, speeding up the checkout and compile times:

<!-- -->

    git clone via ssh:
    bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh)
    git clone via https:
    bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-https.sh)

##### ROOT5 SLC6 release (CMSSW\_7\_1\_X)

**Setting up the environment (once)**

    setenv SCRAM_ARCH slc6_amd64_gcc481
    cmsrel CMSSW_7_1_5 ### must be a 7_1_X release  &gt;= 7_1_5;  (7.0.X and 7.2.X are NOT supported either) 
    cd CMSSW_7_1_5/src 
    cmsenv
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

-   if you get errors related to ZLIB when doing the git clone, try instead doing the git clone before the cmsenv
-   when running with the HybridNew calculator, you can safely ignore the following warning from RooFit
    **`WARNING:Eval -- RooStatsUtils::MakeNuisancePdf - no constraints found on nuisance parameters in the input model`**

**Updating to a tag (both the first time and whenever there are updates)**

    cd HiggsAnalysis/CombinedLimit
    git fetch origin
    git checkout v5.0.4   # try v5.0.1 if any issues occur
    scramv1 b clean; scramv1 b # always make a clean build, as scram doesn't always see updates to src/LinkDef.h

##### SLC5 release (not recommended)

**Setting up the environment (once)**

    setenv SCRAM_ARCH slc5_amd64_gcc472 
    cmsrel CMSSW_6_1_1 ### must be &gt;= 6.1.1, as older versions have bugs (6.2.X and 7.0.X are NOT supported either) 
    cd CMSSW_6_1_1/src 
    cmsenv
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

-   if you get errors related to ZLIB when doing the git clone, try instead doing the git clone before the cmsenv
-   when running with the HybridNew calculator, you can safely ignore the following warning from RooFit
    **`WARNING:Eval -- RooStatsUtils::MakeNuisancePdf - no constraints found on nuisance parameters in the input model`**

**Updating to a tag (both the first time and whenever there are updates)**

    cd HiggsAnalysis/CombinedLimit
    git fetch origin
    git checkout v4.0.1-sl5
    scramv1 b clean; scramv1 b # always make a clean build, as scram doesn't always see updates to src/LinkDef.h

#### %PURPLE%For developers<span class="twiki-macro ENDCOLOR"></span>

We use the *Fork and Pull* model for development: each user creates a copy of the repository on github, commits their requests there and then sends pull requests for the administrators to merge (the admistrators are currently André, Giovanni, Nick and Mingshui )

Prerequisites 1 Register on github, as needed anyway for CMSSW development: <http://cms-sw.github.io/cmssw/faq.html> 1 Register your SSH key on github: <https://help.github.com/articles/generating-ssh-keys> 1 Fork the repository to create your copy of it: <https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/fork> (more documentation at <https://help.github.com/articles/fork-a-repo> ) You will now be able to browse your fork of the repository from <https://github.com/your-github-user-name/HiggsAnalysis-CombinedLimit>

**Setting up the environment (once)**
Get the SSH url of your forked copy from from github. You can find it in the bottom right of the github page for your repository. In my case, it's **`git@github.com:gpetruc/HiggsAnalysis-CombinedLimit.git`**

    setenv SCRAM_ARCH slc6_amd64_gcc481
    cmsrel CMSSW_7_1_5 ### must be a 7_1_X release  &gt;= 7_1_5;  (7.0.X and 7.2.X are NOT supported either) 
    cd CMSSW_7_1_5/src 
    cmsenv
    git clone git@github.com:your-user-name/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit
    git remote add upstream  https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git
    git fetch upstream
    git checkout -b slc6-root5.34.17  upstream/slc6-root5.34.17 

**Recommended way to develop a feature (in a branch)**
<https://blogs.atlassian.com/2013/07/git-upstreams-forks/>

    # get the updates of the master branch of the remote repository
    git fetch upstream

    # branch straight off the upstream master
    git checkout -b feature_name_branch upstream/slc6-root5.34.17 

    # implement the feature
    # commit, etc

    # before publishing:
    # get the updates of the master branch of the remote repository
    git fetch upstream
    # if you're ready to integrate the upstream changes into your repository do
    git rebase upstream/slc6-root5.34.17 
    # fix any conflicts
    git push origin feature_name_branch

And proceed to make a pull request from the branch you created.

**Committing changes to your repository**

    git add ....
    git commit -m "...."
    git push 

**Making a pull request to have your commit/branch integrated upstream**
<https://help.github.com/articles/using-pull-requests> 1 go on your github page: <https://github.com/your-github-user-name/HiggsAnalysis-CombinedLimit> 1 click "pull request" at the right end of the gray bar
&lt;img alt="pull-request.png" height="155" src="%ATTACHURLPATH%/pull-request.png" width="796" /&gt; 1 check that the diff below makes sense to you 1 enter some description of your change and press "create pull request"
&lt;img alt="pull-request-2.png" height="210" src="%ATTACHURLPATH%/pull-request-2.png" width="1004" /&gt;
&lt;img alt="pull-request-3.png" height="368" src="%ATTACHURLPATH%/pull-request-3.png" width="938" /&gt;

### Legacy CVS recipe (NOT FOR DEVELOPMENT)

%BLUE%**Production version**<span class="twiki-macro ENDCOLOR"></span>

    setenv SCRAM_ARCH slc5_amd64_gcc472 
    cmsrel CMSSW_6_1_1 ### must be &gt;= 6.1.1, as older versions have bugs
    cd CMSSW_6_1_1/src 
    cmsenv
    addpkg HiggsAnalysis/CombinedLimit V03-01-12
    scramv1 b

NOTE: when running with the HybridNew calculator, you can safely ignore the following warning from RooFit
**`WARNING:Eval -- RooStatsUtils::MakeNuisancePdf - no constraints found on nuisance parameters in the input model`**

\*Changelog:\*

-   new in **`V03-01-08`** wrt to **`V03-01-06`**: integration of Zγ pdfs and models, and some protection in HybridNew for background models with no parameters
-   new in **`V03-01-06`** wrt to **`V03-01-04`**: fix F-C in HybridNew for one-dimensional models (skipped one version due to a mistake)
-   new in **`V03-01-04`** wrt **`V03-01-03`**: fix a bug introduced in **`V03-01-01`** that prevented running CLs with toys.
-   new in **`V03-01-03`** wrt **`V03-01-02`**: added a feature to create pseudo-Asimov datasets for models with more than one observable (still under testing, will be documented when done)
-   new in **`V03-01-02`** wrt **`V03-01-01`**: fix issues in toy mc generation when using channels a RooDataHist with non-uniform binning as input

<span class="twiki-macro PURPLE"></span> **Old production version** <span class="twiki-macro ENDCOLOR"></span>

```
setenv SCRAM_ARCH slc5_amd64_gcc434 
cmsrel CMSSW_5_2_5   # or any 5.2.X, 5.3.X release
cd CMSSW_5_2_5/src 
cmsenv
addpkg HiggsAnalysis/CombinedLimit V02-07-03
scramv1 b
```

NOTEs:

-   you can safely ignore the message complaining about **`/usr/lib64/libxml2.so.2: no version information available`**
-   %RED%Do **NOT** use this release if you have 2D RooFit histograms with non-uniform binning.<span class="twiki-macro ENDCOLOR"></span>

<span class="twiki-macro RED"></span> **Development version** <span class="twiki-macro ENDCOLOR"></span>

```
setenv SCRAM_ARCH slc5_amd64_gcc472 
cmsrel CMSSW_6_1_1
cd CMSSW_6_1_1/src 
cmsenv
cvs co HiggsAnalysis/CombinedLimit
scramv1 b
```

NOTE: you can safely ignore the warning from RooFit
**`WARNING:Eval -- RooStatsUtils::MakeNuisancePdf - no constraints found on nuisance parameters in the input model`**

Note: If you're using the **`tcsh`** shell instead of the **`bash`** one, you should also execute the command **`rehash`** after you have compiled the program.

The package contains one binary executable **`combine`** which will take as input one datacard describing your model and compute a cross section limit or the significance of an observation.

The package can also be combined against ROOT in a stand-alone mode (provided that Boost is available), although this is not officially supported. This can be achieved as follows:

```
cvs co -d ./ HiggsAnalysis/CombinedLimit
cd CombinedLimit
make
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/lib
```

If your Boost includes are not in the default directory, you can edit the first line in the makefile so to specify the correct path.

@subpage MLFits