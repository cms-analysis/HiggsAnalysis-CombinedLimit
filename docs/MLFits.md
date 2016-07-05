Maximum likelihood fits and diagnostics {#MLFits}
=================================================

If you have already found the higgs boson but it's an exotic one, instead of computing a limit or significance you might want to extract it's cross section by performing a maximum-likelihood fit. Or, more seriously, you might want to use this same package to extract the cross section of some other process (e.g. the di-boson production). Or you might want to know how well the data compares to you model, e.g. how strongly it constraints your other nuisance parameters, what's their correlation, ...

Running the MaxLikelihoodFit tool as

    combine -M MaxLikelihoodFit datacard.txt 

The program will print out the result of the two fits performed with signal strength set to zero and with floating signal strength, and the output root tree will contain the best fit value for the signal strength, and it's uncertainty. You will also get a `mlfit.root` file containing the following objects:

|                        |                                                                                                                                       |
|------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| **`nuisances_prefit`** | `RooArgSet` containing the pre-fit values of the nuisance parameters, and their uncertainties from the external constraint terms only |
| **`fit_b`**            | `RooFitResult` object containing the outcome of the fit of the data with signal strength set to zero                                  |
| **`fit_s`**            | `RooFitResult` object containing the outcome of the fit of the data with floating signal strength                                     |
| **`covariance_fit_s`** | `TH2D` Covariance matrix of the parameters in the fit with signal strength set to zero                                                |
| **`covariance_fit_b`** | `TH2D` Covariance matrix of the parameters in the fit with floating signal strength                                                   |

It is possible to compare pre-fit and post-fit nuisance parameters with the script [diffNuisances.py](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/test/diffNuisances.py). Taking as input a `mlfit.root` file, the script will by default print out the parameters which have changed significantly w.r.t. their initial estimate. For each of those parameters, it will print out the shift in value and the post-fit uncertainty, both normalized to the input values, and the linear correlation between the parameter and the signal strength.

    python diffNuisances.py  mlfit.root 

The script has several options to toggle the thresholds used to decide if a parameter has changed significantly, to get the printout of the absolute value of the nuisance parameters, and to get the output in another format for easy cut-n-paste (supported formats are `html`, `latex`, `twiki`). In addition, this script has the option (`-g outputfile.root`) to produce plots of the fitted values of the nuisance parameters and their post-fit uncertainties as well as a plot showing directly a comparison of the post-fit to pre-fit nuisance uncertainties.

**Background normalizations**
For a certain class of models, like those made from datacards for shape-based analysis, the tool can also compute and save to the output root file the best fit yields of all background processes. If this feature is turned on with the option **`--saveNormalizations`**, the output root file will contain also two RooArgSet `norm_fit_s`, `norm_fit_b` objects each containing one RooConstVar for each channel `xxx` and process `yyy` with name `xn_exp_bin<xxx>_proc_<xxx>` and value equal to the best fit yield. In the future, also the uncertainties on those yields will be reported.
Note that this procedure works only for "extended likelihoods" like the ones used in shape-based analysis, not for the cut-and-count datacards. You can however convert a cut-and-count datacard in an equivalent shape-based one by adding a line `shapes * * FAKE` in the datacard after the `imax`, `jmax`, `kmax`.

The sample pyroot macro [mlfitNormsToText.py](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/test/mlfitNormsToText.py) can be used to convert the root file into a text table with four columns: channel, process, yield from the signal+background fit and yield from the background-only fit.

**Plots**
MaxLikelihoodFit can also produce pre- and post-fit plots of data and input model in the same directory as `mlfit.root`. To get them, you have to specify the option **`--plots`**, and then *optionally specify* what are the names of the signal and background pdfs, e.g. **`--signalPdfNames='ggH*,vbfH*'`** and `=--backgroundPdfNames='*DY*,*WW*,*Top*'` (by default, the definitions of signal and background are taken from the datacard). In this case you might want to use the option **`--out`** to specify the output directory where to write files. <span class="twiki-macro RED"></span> Currently this is very very crude, **and in some cases (eg binned pdfs with different bin widths) will give incorrect results so their interpretation should be taken with extreme caution**. If people find it useful they could consider contributing to its development.

An alternative is to use the options `--saveShapes` (`--saveWithUncertainties` for adding post-fit uncertainties). The result will be additional folders in `mlfit.root` for each category, with pre and post-fit distributions of the signals and backgrounds (as TH1). These distributions and normalisations are guaranteed to give the correct interpretation however, the data is not included. This can be taken of course from the input workspaces (adding a histogram for the data is currently in development). It should be noted that for shape datacards whose inputs are TH1, these histograms will have the bin number as the x-axis and the content of each bin will be a number of events. Instead for datacards whose inputs are RooAbsPdf, the x-axis will correspond to the observable and the bin content will be the PDF density (i.e the number of events in a given bin, i, can be obtained from h-&gt;GetBinContent(i)\*h-&gt;GetBinWidth(i) ).

##### Fit options

-   If you need only the signal+background fit, you can run with **`--justFit`**. This can be useful if the background-only fit is not interesting or not converging (e.g. if the significance of your signal is very very large)
-   You can use **`--rMin`** and **`--rMax`** to set the range of the parameter; a range that is not too large compared to the uncertainties you expect from the fit usually gives more stable and accurate results.
-   By default, the uncertainties are computed using MINOS; if this fails to converge, you can try running with **`--robustFit=1`** that will do a slower but more robust likelihood scan; this can be further controlled by the parameter **`--stepSize`** (the default is 0.1, and is relative to the range of the parameter)
-   Different configurations of the fitter can be tested to debug convergence problems:
    -   **`--minimizerAlgo`** can be set to `Minuit2` (default) or `Minuit`
    -   **`--minimizerTolerance`** can be varied between 0.0001 and 1.0 or more (the default is 0.01); smaller tolerances can give a more accurate result and avoid false minima, while larger tolerances can allow a complex fit to converge to a satisfactory result rather than giving up.
    -   **`--minimizerStrategy`** can be set to 1 (default), 0 or 2; larger values correspond to using more derivatives in the fit, which is slower but can give better or more stable results for smooth pdfs (and can result in failures for less smooth pdfs)

##### Toy-by-toy diagnostics

MaxLikelihoodFit can also be used to diagnose the fitting procedure in toy experiments to identify potentially problematic nuisance parameters when running the full limits/p-values. This can be done by adding the option -t Ntoys. The output file will contain two trees, one for the background only fit and one for the signal + background, which contain the value of the constraint fitted result in each toy. It is recommended to use the following options when investigating toys to reduce the running time: **`--toysFrequentist`** **`--noErrors`** **`--minos none`**

The results can be plotted using plotParametersFromToys.C.

    plotParametersFomToys("mlfittoys.root","mlfitdata.root","workspace.root","mu<0")

The first argument is the name of the output file from running with toys, and the second and third (optional) arguments are the name of the file containing the result from a fit to the data and the workspace. The fourth argument can be used to specify a cut string applied to one of the branches in the tree which can be used to correlate strange behaviour with specific conditions.

##### Calculate expected precision on signal strength

In order to estimate the error on the signal strength using MC information only (Asimov dataset):

    combine -M MaxLikelihoodFit --expectSignal=1 -t -1 datacard.txt

-   **`expectSignal`** is the simulated signal strength
-   **`-t -1`** is to use the Asimov dataset
-   If interested only in the the outcome of the s+b fit, then add an extra **`--justFit`** which will skip the background-only fit and the rest of the diagnostics.

#### Computing the expected limit (bayesian or likelihood methods)

For bayesian or likelihood-based limits, the expected limit is computed by generating many toy mc observations and compute the limit for each of them. This can be done passing the option **`-t`** . E.g. to run 100 toys with the ProfileLikelihood method, just do

    combine -M ProfileLikelihood datacard.txt -t 100

The program will print out the mean and median limit, and the 68% and 95% quantiles of the distributions of the limits. Also, the output root tree will contain one entry per toy.

While for a simple method and the ProfileLikelihood method it can be feasible to run all the toys in a single job, for more heavy methods you'll probably want to split this in multiple jobs. To do this, just run `combine` multiple times specifying a smaller number of toys (can be as low as `1`) each time using a different seed to initialize the random number generator (option `=-s==; if you set it to -1, the starting seed will be initialized randomly at the beginning of the job), then merge the resulting trees with =hadd` and look at the distribution in the merged file.

Macros to turn the rootfiles into bands will be committed to the cvs soon.

#### Computing the expected limit (hybrid or frequentist methods)

There are two possible approaches for computing expected limits with HybridNew, depending on the the complexity of the model.

-   Run interactively 5 times to compute the median expected and the +/- 1 and +/- 2 sigma lines
-   Produce once a grid of test statistics distributions at various values of the signal strength, and use it to compute the observed and expected limit and bands

The first approach is good for simple models, the second is good for complex models since the grid of points can be distributed across any number of jobs.

For the first approach, just run interactively HybridNew with the same options as per the observed limit but adding a **`--expectedFromGrid=<quantile>`** where the quantile is 0.5 for the median, 0.84 for the +1 sigma band, 0.16 for the -1 sigma band, 0.975 for the +2 sigma band, 0.025 for the -2 sigma band.

For the second approach, read under [HybridNew and grids](%SCRIPTURL{"view"}%auth/CMS/SWGuideHiggsAnalysisCombinedLimit#HybridNew_algorithm_usage_for_co).

#### Computing the expected significance (from asymptotic)

The expected significance can be computed either from an Asimov dataset of signal+background or from toys of (signal+background). For both approaches, there are two alternative definitions:

-   a positeriori expected: depends on the observed dataset.
-   a priori expected: does not depend on the observed dataset, and so is a good metric for optimizing an analysis when still blinded.

The a-posteriori expected significance from the Asimov dataset, which is what was used for the Higgs search results in 2012, is

    combine -M ProfileLikelihood --significance datacard.txt -t -1 --expectSignal=1 --toysFreq

the a-priori expectation can be computed in a very similar way by omitting the **`--toysFreq`** argument

    combine -M ProfileLikelihood --significance datacard.txt -t -1 --expectSignal=1

The a-posteriori and a-priori expected significances with the toys can be computed in the similar way but replacing **`-t -1`** with **`-t N`** (N being the number of toys to run). As for other toy-based computations, the task can be split in multiple jobs run with different random number generator seeds (option **`-s`**)

The output format is the same as for expected limits: the variable **`limit`** in the tree will be filled with the significance (or with the p-value if you put also the option **`--pvalue`**)

#### Computing the expected significance (toys)

Computation of expected significance with toys is a two step procedure: first you need to run one or more jobs to construct the expected distribution of the test statistics, e.g.

    combine -M HybridNew --frequentist datacard.txt --significance  [ --pvalue ] --saveToys --fullBToys --saveHybridResult   -T toys -i iterations -s seed 

then you can merge all those results into a single root file with **`hadd`** and compute the median expected significance as

    combine -M HybridNew --frequentist datacard.txt --significance  [ --pvalue ] --readHybridResult --toysFile=input.root --expectedFromGrid=0.5

to get edges the ±1σ bands, use 0.16 and 0.84 instead of 0.5, and so on...

Note that you need a total number of background toys large enough to compute the value of the significance, but you need less signal toys (especially if you only need the median). For large significance, you can then run most of the toys without the **`--fullBToys`** limit (about a factor 2 faster), and only a smaller part with that option turned on.
