# Useful links and further reading

### Tutorials and reading material

There are several tutorials which have been run over the last few years with instructions and examples for running the combine tool.

* Tutorial Sessions
   * [1st tutorial 17th Nov 2015](https://indico.cern.ch/event/456547/).
   * [2nd tutorial 30th Nov 2016](https://indico.cern.ch/event/577649/#b-229590-higgs-combine-tool-mi).
   * [3rd tutorial 29th Nov 2017](https://indico.cern.ch/event/677948/#day-2017-11-29) 
   * [4th tutorial 31st Oct 2018](https://indico.cern.ch/event/747340/overview) - Latest for `81x-root606` branch.
   * [5th tutorial 2nd-4th Dec 2019](https://indico.cern.ch/event/859454/overview)
   * [6th tutorial 14th-16th Dec 2020](https://indico.cern.ch/event/976099/overview) - Latest for `102x` branch


* Worked examples from Higgs analyses using combine  
   * [The CMS DAS at CERN 2014](https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideCMSDataAnalysisSchool2014HiggsCombPropertiesExercise)
   * [The CMS DAS at DESY 2018](https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchoolHamburg2018LongStatisticsExercise)



* [Conventions to be used when preparing inputs for Higgs combinations](https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationConventions)

* [CMS AN-2011/298](http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS AN-2011/298) Procedure for the LHC Higgs boson search combination in summer 2011. This describes in more detail some of the methods used in Combine.

### Citations

There is no document currently which can be cited for using the combine tool, however you can use the following publications for the procedures we use,

* [Summer 2011 public ATLAS-CMS note](https://cds.cern.ch/record/1379837) for any Frequentist limit setting procedures with toys or Bayesian limits, constructing likelihoods, descriptions of nuisance parameter options (like log-normals (`lnN`) or gamma (`gmN`), and for definitions of test-statistics.

* [CCGV paper](https://arxiv.org/abs/1007.1727) if you use any of the **asymptotic** (eg with `-M AsymptoticLimits` or `-M Significance` approximations for limits/p-values.

* If you use the Barlow-Beeston approach to MC stat (bin-by-bin) uncertainties, please cite their paper [Barlow-Beeston](https://www.sciencedirect.com/science/article/pii/001046559390005W?via%3Dihub). You should also cite [this note](https://arxiv.org/pdf/1103.0354.pdf) if you use the `autoMCStats` directive to produce a single parameter per bin.

* If you use `shape` uncertainties for template (`TH1` or `RooDataHist`) based datacards, you can cite [this note](https://arxiv.org/pdf/1103.0354.pdf) from J. Conway.

* If you are extracting uncertainties from LH scans - i.e using $-2\Delta Log{L}=1$ etc for the 1$\sigma$ intervals, you can cite either the [ATLAS+CMS](https://arxiv.org/abs/1606.02266) or [CMS](https://link.springer.com/article/10.1140/epjc/s10052-015-3351-7) Higgs paper.

* There is also a long list of citation recommendations from the [CMS Statistics Committee](https://twiki.cern.ch/twiki/bin/view/CMS/StatisticsReferences) pages.

### Combine based packages

* [SWGuideHiggs2TauLimits](https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggs2TauLimits)

* [ATGCRooStats](https://twiki.cern.ch/twiki/bin/view/CMS/ATGCRooStats)

* [CombineHarvester](http://cms-analysis.github.io/CombineHarvester/)

### Contacts

* **Hypernews forum**: hn-cms-higgs-combination [https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination.html](https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination.html)

### CMS Statistics Committee

* You can find much more statistics theory and reccomendations on various statistical procedures in the [CMS Statistics Committee Twiki Pages](https://twiki.cern.ch/twiki/bin/viewauth/CMS/StatisticsCommittee#Recommendations_from_the_Committ)

# FAQ

* _Why does combine have trouble with bins that have zero expected contents?_
    * If you're computing only upper limits, and your zero-prediction bins are all empty in data, then you can just set the background to a very small value instead of zero as anyway the computation is regular for background going to zero (e.g. a counting experiment with $B\leq1$ will have essentially the same expected limit and observed limit as one with $B=0$). If you're computing anything else, e.g. p-values, or if your zero-prediction bins are not empty in data, you're out of luck, and you should find a way to get a reasonable background prediction there (and set an uncertainty on it, as per the point above)
* _How can an uncertainty be added to a zero quantity?_
    * You can put an uncertainty even on a zero event yield if you use a gamma distribution. That's in fact the more proper way of doing it if the prediction of zero comes from the limited size of your MC or data sample used to compute it.
* _Why does changing the observation in data affect my expected limit?_
    * The expected limit (if using either the default behaviour of `-M AsymptoticLimits` or using the `LHC-limits` style limit setting with toys) uses the _**post-fit**_ expectation of the background model to generate toys. This means that first the model is fit to the _**observed data**_ before toy generation. See the sections on [blind limits](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/#blind-limits) and [toy generation](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/runningthetool/#toy-data-generation) to avoid this behavior. 
* _How can I deal with an interference term which involves a negative contribution?_
    * You will need to set up a specific PhysicsModel to deal with this, however you can [see this section](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part2/physicsmodels/#interference) to implement such a model which can incorperate a negative contribution to the physics process
* _How does combine work?_
    * That is not a question which can be answered without someone's head exploding so please try to formulate something specific.
* _What does fit status XYZ mean?_ 
    * Combine reports the fit status in some routines (for example in the `FitDiagnostics` method). These are typically the status of the last call from Minuit. For details on the meanings of these status codes see the [Minuit2Minimizer](https://root.cern.ch/root/html/ROOT__Minuit2__Minuit2Minimizer.html) documentation page.
* _Why does my fit not converge?_ 
    * There are several reasons why some fits may not converge. Often some indication can be obtained from the `RooFitResult` or status which you will see information from when using the `--verbose X` (with $X>2$) option. Sometimes however, it can be that the likelihood for your data is very unusual. You can get a rough idea about what the likelihood looks like as a function of your parameters (POIs and nuisances) using `combineTool.py -M FastScan -w myworkspace.root` (use --help for options). 
    * We have  seen often that fits in combine using `RooCBShape` as a parametric function will fail. This is related to an optimisation that fails. You can try to fix the problem as described in this issue: https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/issues/347 (i.e add the option `--X-rtd ADDNLL_CBNLL=0`).
* _Why does the fit/fits take so long?_ 
    * The minimisation routines are common to many methods in combine. You can tune the fitting using the generic optimisation command line options described [here](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/runningthetool/#generic-minimizer-options). For example, setting the default minimizer strategy to 0 can greatly improve the speed since this avoids running Hesse. In calculations such as `AsymptoticLimits`, Hesse is not needed and hence this can be done, however, for `FitDiagnostics` the uncertainties and correlations are part of the output so using strategy 0 may not be particularly accurate. 
* _Why are the results for my counting experiment so slow or unstable?_ 
    * There is a known issue with counting experiments with ***large*** numbers of events which will cause unstable fits or even the fit to fail. You can avoid this by creating a "fake" shape datacard (see [this section](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part2/settinguptheanalysis/#combination-of-multiple-datacards) from the setting up the datacards page). The simplest way to do this is to run `combineCards.py -S mycountingcard.txt > myshapecard.txt`. You may still find that your parameter uncertainties are not correct when you have large numbers of events. This can be often fixed using the `--robustHesse` option. An example of this issue is detailed [here](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/issues/498). 
* _Why do some of my nuisance parameters have uncertainties &gt; 1?_
    * When running `-M FitDiagnostics` you may find that the post-fit uncertainties of the nuisances are $> 1$ (or larger than their pre-fit values). If this is the case, you should first check if the same is true when adding the option `--minos all` which will invoke minos to scan the likelihood as a function of these parameters to determine the crossing at $-2\times\Delta\log\mathcal{L}=1$ rather than relying on the estimate from Hesse. However, this is not guaranteed to succeed, in which case you can scan the likelihood yourself using `MultiDimFit` (see [here](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/#likelihood-fits-and-scans) ) and specifying the option `--poi X` where `X` is your nuisance parameter. 
* _How can I avoid using the data?_ 
    * For almost all methods, you can use toy data (or an Asimov dataset) in place of the real data for your results to be blind. You should be careful however as in some methods, such as `-M AsymptoticLimits` or `-M HybridNew --LHCmode LHC-limits` or any other method using the option `--toysFrequentist`, the data will be used to determine the most likely nuisance parameter values (to determine the so-called a-posteriori expectation). See the section on [toy data generation](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/runningthetool/#toy-data-generation) for details on this. 
* _What if my nuisance parameters have correlations which are not 0 or 1?_
    * Combine is designed under the assumption that each *source* of nuisance parameter is uncorrelated with the other sources. If you have a case where some pair (or set) of nuisances have some known correlation structure, you can compute the eigenvectors of their correlation matrix and provide these *diagonalised* nuisances to combine. You can also model *partial correlations*, between different channels or data taking periods, of a given nuisance parameter using the `combineTool` as described in [this page](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/issues/503). 
* _My nuisances are (artificially) constrained and/or the impact plot show some strange behaviour, especially after including MC statistical uncertainties. What can I do?_
    * Depending on the details of the analysis, several solutions can be adopted to mitigate these effects. We advise to run the validation tools at first, to identify possible redundant shape uncertainties that can be safely eliminated or replaced with lnN ones. Any remaining artificial constrain should be studies. Possible mitigating strategies can be to (a) smooth the templates or (b) adopt some rebinning in order to reduce statistical fluctuations in the templates. A description of possible strategies and effects can be found in [this talk by Margaret Eminizer](https://indico.cern.ch/event/788727/contributions/3401374/attachments/1831680/2999825/higgs_combine_4_17_2019_fitting_details.pdf)
* _What do CLs, CLs+b and CLb in the code mean?_
    * The names CLs+b and CLb are rather outdated and should instead be referred to as p-values - $p_{\mu}$ and $1-p_{b}$, respectively. We use the CLs (which itself is not a p-value) criterion often in High energy physics as it is designed to avoid excluding a signal model when the sensitivity is low (and protects against excluding due to underfluctuations in the data). Typically, when excluding a signal model the p-value $p_{\mu}$ often refers to the p-value under the signal+background hypothesis, assuming a particular value of the signal stregth ($\mu$) while $p_{b}$ is the p-value under the background only hypothesis. You can find more details and definitions of the CLs criterion and $p_{\mu}$ and $p_{b}$ in section 39.4.2.4 the [2016 PDG review](http://pdg.lbl.gov/2016/reviews/rpp2016-rev-statistics.pdf).       


