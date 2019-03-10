# Setting up the environment and installation

The instructions below are for installation within a CMSSW environment

## For end users that don't need to commit or do any development

You can find the latest releases on github under [https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/releases)

### ROOT6 SLC6 release `CMSSW_8_1_X` - recommended version

Setting up the environment (once):

```sh
export SCRAM_ARCH=slc6_64_gcc530
cmsrel CMSSW_8_1_0
cd CMSSW_8_1_0/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
```
Update to a reccomended tag - currently the reccomended tag is **v7.0.12**:

```sh
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v7.0.12
scramv1 b clean; scramv1 b # always make a clean build
```

You can generate a diff of any two tags (eg for `v7.0.8` and `v7.0.6`) by using following the url:

[https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/compare/v7.0.6...v7.0.7](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/compare/v7.0.6...v7.0.7)

Replace the tag names in the url to any tags you which to compare.

## For developers

We use the _Fork and Pull_ model for development: each user creates a copy of the repository on github, commits their requests there and then sends pull requests for the administrators to merge.

_Prerequisites_

1. Register on github, as needed anyway for CMSSW development: [http://cms-sw.github.io/cmssw/faq.html](http://cms-sw.github.io/cmssw/faq.html)

2. Register your SSH key on github: [https://help.github.com/articles/generating-ssh-keys](https://help.github.com/articles/generating-ssh-keys) 1 Fork the repository to create your copy of it: [https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/fork](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/fork) (more documentation at [https://help.github.com/articles/fork-a-repo](https://help.github.com/articles/fork-a-repo) )

You will now be able to browse your fork of the repository from [https://github.com/your-github-user-name/HiggsAnalysis-CombinedLimit](https://github.com/your-github-user-name/HiggsAnalysis-CombinedLimit)

### Recommended way to develop a feature (in a branch)

```sh
# get the updates of the master branch of the remote repository
git fetch upstream

# branch straight off the upstream master
git checkout -b feature_name_branch upstream/81x-root606

# implement the feature
# commit, etc

# before publishing:
# get the updates of the master branch of the remote repository
git fetch upstream

# if you're ready to integrate the upstream changes into your repository do
git rebase upstream/81x-root606

# fix any conflicts
git push origin feature_name_branch
```

And proceed to make a pull request from the branch you created.

#### Committing changes to your repository

```sh
git add ....
git commit -m "...."
git push
```

You can now make a pull request to the repository.

## Combine Tool

An additional tool for submitting combine jobs to batch/crab, developed originally for HiggsToTauTau. Since the repository contains a certain amount of analysis-specific code, the following scripts can be used to clone it with a sparse checkout for just the core [`CombineHarvester/CombineTools`](https://github.com/cms-analysis/CombineHarvester/blob/master/CombineTools/) subpackage, speeding up the checkout and compile times:

git clone via ssh:

```sh
bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh)
```

git clone via https:

```sh
bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-https.sh)
```

make sure to run `scram`  to compile the `CombineTools` package.

See the [`CombineHarvester`](http://cms-analysis.github.io/CombineHarvester/) documentation pages for more details on using this tool and additional features available in the full package.

