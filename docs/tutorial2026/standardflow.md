# Main Features of Combine
This exercise is designed to recreate the main workflow needed to perform a statistical analysis with Combine. It will start assuming you already prepared your inputs (**shapes, yields, and systematic uncertainties**) and will proceed step by step to perform validation test of your setup and produce some standard results. For more detailed procedure you can always find detailed informations in the Combine **manual** and in the **Long exercise tutorial**. 

As for the Long exercise, we will work with a simplified version of a real analysis, that nonetheless will have many features of the full analysis. The analysis is a search for an additional heavy neutral Higgs boson decaying to tau lepton pairs. Such a signature is predicted in many extensions of the standard model, in particular the minimal supersymmetric standard model (MSSM). You can read about the analysis in the paper [here](https://arxiv.org/pdf/1803.06553.pdf). The statistical inference makes use of a variable called the total transverse mass ($M_{\mathrm{T}}^{\mathrm{tot}}$) that provides good discrimination between the resonant high-mass signal and the main backgrounds, which have a falling distribution in this high-mass region. The events selected in the analysis are split into a several categories which target the main di-tau final states as well as the two main production modes: gluon-fusion (ggH) and b-jet associated production (bbH). One example is given below for the fully-hadronic final state in the b-tag category which targets the bbH signal:

![](images/CMS-DAS.003.jpeg)

## Background
You can find a presentation with some more background on likelihoods and extracting confidence intervals [here](https://indico.cern.ch/event/976099/contributions/4138517/). A presentation that discusses limit setting in more detail can be found [here](https://indico.cern.ch/event/976099/contributions/4138520/).
If you are not yet familiar with these concepts, or would like to refresh your memory, we recommend that you have a look at these presentations before you start with the exercise.

## Getting started
To get started, you should have a working setup of `Combine`, please follow the instructions from the [home page](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#within-cmssw-recommended-for-cms-users). Make sure to use the latest recommended release.

Now we will move to the working directory for this tutorial, which contains all the inputs needed to run the exercises below:
```shell
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/docs/tutorial2026/
```


