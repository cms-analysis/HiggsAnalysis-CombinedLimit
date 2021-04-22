# How to run the tool

The executable **`combine`** provided by the package allows to use the Higgs Combination Tool indicating by command line which is the method to use for limit combination and which are user's preferences to run it. To see the entire list of all available options ask for the help:

```sh
combine --help
```

The option `-M` allows to chose the method used. There are several groups of statistical methods:

-   **Asymptotic** likelihood methods:
    -   `AsymptoticLimits`: limits calculated according to the asymptotic formulas in [arxiv:1007.1727](http://arxiv.org/abs/1007.1727)
    -   `Significance`: simple profile likelihood approximation, for calculating significances.
-   **Bayesian** methods:
    -   `BayesianSimple`: performing a classical numerical integration (for simple models only)
    -   `MarkovChainMC`: performing Markov Chain integration, for arbitrarily complex models.
-   **Frequentist** or hybrid bayesian-frequentist methods:
    -   `HybridNew`: compute modified frequentist limits according to several possible prescriptions
-   **Fitting**
    - `FitDiagnostics`: performs maximum likelihood fits to extract the signal yield and provide diagnostic tools such as pre and post-fit models and correlations
    - `MultiDimFit`: perform maximum likelihood fits in multiple parameters and likelihood scans
-   **Miscellaneous** other modules that don't compute limits but use the same framework:
    - `GoodnessOfFit`: perform a goodness of fit test for models including shape information using several GOF estimators
    - `ChannelConsistencyCheck`: check how consistent are the individual channels of a combination are
    - `GenerateOnly`: generate random or asimov toy datasets for use as input to other methods

The command help is organized into five parts:

-   *Main options* section indicates how to pass the datacard as input to the tool (`-d datacardName`) and how to choose the statistical method (`-M MethodName`) to compute a limit and level of verbosity for output `-v`
-   *Common statistics options* include options common to different statistical methods such as `--cl` to specify the CL (default is 0.95) or `-t` to give the number of toy MC extractions required.
-   *Common input-output options*. Is it possible to specify hypothesis point under analysis using `-m` or include specific string in output filename `--name`.
-   *Common miscellaneous options*.
-   Method specific options sections are dedicated to each method. By providing the Method name with the `-M` option, only the options for that specific method are shown in addition to the common options

Those options reported above are just a sample of all available.The command `--help` provides documentation of all of them.

## Common command line options

There are a number of useful command line options which can be used to alter the model (or parameters of the model) at run These are the most commonly used, generic options,

-   `-H`: run first another faster algorithm (e.g. the ProfileLikelihood described below) to get a hint of the limit, allowing the real algorithm to converge more quickly. We **strongly recommend** to use this option when using MarkovChainMC, HybridNew or FeldmanCousins calculators, unless you know in which range your limit lies and you set it manually (the default is `[0, 20]`)

-   `--rMax`, `--rMin`: manually restrict the range of signal strengths to consider. For Bayesian limits with MCMC, `rMax` a rule of thumb is that rMax should be 3-5 times the limit (a too small value of `rMax` will bias your limit towards low values, since you are restricting the integration range, while a too large value will bias you to higher limits)

-   `--setParameters name=value[,name2=value2,...]` sets the starting values of the parameters, useful e.g. when generating toy MC or when also setting the parameters as fixed. This option supports the use of regexp via by replacing `name` with `rgx{some regular expression}`.

-   `--setParameterRanges name=min,max[:name2=min2,max2:...]` sets the ranges of the parameters (useful e.g. for scanning in MultiDimFit, or for Bayesian integration). This option supports the use of regexp via by replacing `name` with `rgx{some regular expression}`.

-   `--redefineSignalPOIs name[,name2,...]` redefines the set of parameters of interest.
    -   if the parameters where constant in the input workspace, they are re-defined to be floating.
    -   nuisances promoted to parameters of interest are removed from the list of nuisances, and thus they are not randomized in methods that randomize nuisances (e.g. HybridNew in non-frequentist mode, or BayesianToyMC, or in toy generation with `-t` but without `--toysFreq`). This doesn't have any impact on algorithms that don't randomize nuisances (e.g. fits, AsymptoticLimits, or HybridNew in fequentist mode) or on algorithms that treat all parameters in the same way (e.g. MarkovChainMC).
    -   Note that constraint terms for the nuisances are **dropped** after promotion to a POI using `--redefineSignalPOI`. To produce a likelihood scan for a nuisance parameter, using MultiDimFit with **`--algo grid`**, you should instead use the `--parameters (-P)` option which will not cause the loss of the constraint term when scanning.
    -   parameters of interest of the input workspace that are not selected by this command become unconstrained nuisance parameters, but they are not added to the list of nuisances so they will not be randomized (see above).


-   `--freezeParameters name1[,name2,...]` Will freeze the parameters with the given names to their set values. This option supports the use of regexps via by replacing `name` with `rgx{some regular expression}` for matching to *constrained nuisance parameters* or `var{some regular expression}` for matching to *any* parameter. For example `--freezeParameters rgx{CMS_scale_j.*}` will freeze all constrained nuisance parameters with the prefix `CMS_scale_j`, while `--freezeParameters var{.*rate_scale}` will freeze any parameter (constrained nuisance or otherwise) with the suffix `rate_scale`.
    - use the option `--freezeParameters allConstrainedNuisances` to freeze all nuisance parameters that have a constraint term (i.e not `flatParams` or `rateParams` or other freely floating parameters).
    - similarly the option `--floatParameters` sets the parameter floating.
    - groups of nuisances (constrained or otherwise), as defined in the datacard, can be frozen using `--freezeNuisanceGroups`. You can also specify to freeze nuisances which are *not* contained in a particular group using a **^** before the group name (`--freezeNuisanceGroups=^group_name` will freeze everything except nuisance parameters in the group "group_name".)
    - all *constrained* nuisance parameters (not `flatParam` or `rateParam`) can be set floating using `--floatAllNuisances`.

!!! warning
    Note that the floating/freezing options have a priority ordering from lowest to highest as `floatParameters < freezeParameters < freezeNuisanceGroups < floatAllNuisances`. Options with higher priority will override those with lower priority.

-   `--trackParameters name1[,name2,...]` will add a branch to the output tree for each of the named parameters. This option supports the use of regexp via by replacing `name` with `rgx{some regular expression}`

    - the name of the branch will be **trackedParam_*name***.
    - the exact behaviour depends on the method. For example, when using `MultiDimFit` with the `--algo scan`, the value of the parameter at each point in the scan will be saved while for `FitDiagnostics`, only the value at the end of the method will be saved.

-   `--trackErrors name1[,name2,...]` will add a branch to the output tree for the error of each of the named parameters. This option supports the use of regexp via by replacing `name` with `rgx{some regular expression}`

    - the name of the branch will be **trackedError_*name***.
    - the behaviour is the same as `--trackParameters` above.

#### Generic Minimizer Options

Combine uses its own minimizer class which is used to steer Minuit (via RooMinimizer) named the `CascadeMinimizer`. This allows for sequential minimization which can help in case a particular setting/algo fails. Also, the `CascadeMinimizer` knows about extra features of Combine such as *discrete* nuisance parameters.

All of the fits which are performed in several of the methods available use this minimizer. This means that the fits can be tuned using these common options,

* `--cminPoiOnlyFit`: First, perform a fit floating *only* the parameters of interest. This can be useful to find, roughly, where the global minimum is.
* `--cminPreScan`: Do a scan before first minimization
* `--cminPreFit` arg: If set to a value N > 0, the minimizer will perform a pre-fit with strategy (N-1) with frozen nuisance parameters.
     * `--cminApproxPreFitTolerance arg`: If non-zero, do first a pre-fit with this tolerance (or 10 times the final tolerance, whichever is largest)
     * `--cminApproxPreFitStrategy arg`:   Strategy to use in the pre-fit. The default is strategy 0.
* `--cminDefaultMinimizerType arg`: Set the default minimizer Type. Default is Minuit2.
* `--cminDefaultMinimizerAlgo arg`: Set the default minimizer Algo. The default is Migrad
* `--cminDefaultMinimizerTolerance arg`: Set the default minimizer Tolerance, the default is 0.1
* `--cminDefaultMinimizerStrategy arg`: Set the default minimizer Strategy between 0 (speed), 1 (balance - *default*), 2 (robustness). The [Minuit documentation](http://www.fresco.org.uk/minuit/cern/node6.html) for this is pretty sparse but in general, 0 means evaluate the function less often, while 2 will waste function calls to get precise answers. An important note is that Hesse (error/correlation estimation) will be run *only* if the strategy is 1 or 2.
* `--cminFallbackAlgo arg`: Provides a list of fallback algorithms if the default minimizer fails. You can provide multiple ones using the syntax is `Type[,algo],strategy[:tolerance]`: eg `--cminFallbackAlgo Minuit2,Simplex,0:0.1` will fall back to the simplex algo of Minuit2 with strategy 0 and a tolerance 0.1, while `--cminFallbackAlgo Minuit2,1` will use the default algo (migrad) of Minuit2 with strategy 1.
* `--cminSetZeroPoint (0/1)`: Set the reference of the NLL to 0 when minimizing, this can help faster convergence to the minimum if the NLL itself is large. The default is true (1), set to 0 to turn off.

The allowed combinations of minimizer types and minimizer algos are as follows

| **Minimizer Type** | **Minimizer Algo** |
|--------------------|--------------------|
|`Minuit`	      | `Migrad`, `Simplex`, `Combined`, `Scan` |
|`Minuit2` 	      | `Migrad`, `Simplex`, `Combined`, `Scan` |
|`GSLMultiMin`       | `ConjugateFR`, `ConjugatePR`, `BFGS`, `BFGS2`, `SteepestDescent`|

More of these options can be found in the **Cascade Minimizer options** section when running `--help`.


#### Output from combine

Most methods will print the results of the computation to the screen, however, in addition, combine will also produce a root file containing a tree called **limit** with these results. The name of this file will be of the format,

    higgsCombineTest.MethodName.mH$MASS.[word$WORD].root

where **$WORD** is any user defined keyword from the datacard which has been set to a particular value.

A few command line options of combine can be used to control this output:

-   The option `-n` allows you to specify part of the name of the rootfile. e.g. if you do `-n HWW` the roofile will be called `higgsCombineHWW....` instead of `higgsCombineTest`
-   The option `-m` allows you to specify the higgs boson mass, which gets written in the filename and also in the tree (this simplifies the bookeeping because you can merge together multiple trees corresponding to different higgs masses using `hadd` and then use the tree to plot the value of the limit vs mass) (default is m=120)
-   The option `-s` allows to specify the seed (eg `-s 12345`) used in toy generation. If this option is given, the name of the file will be extended by this seed, eg `higgsCombineTest.AsymptoticLimits.mH120.12345.root`
-   The option `--keyword-value` allows you to specify the value of a keyword in the datacard such that **$WORD** (in the datacard) will be given the value of **VALUE** in the command `--keyword-value WORD=VALUE`, eg  `higgsCombineTest.AsymptoticLimits.mH120.WORDVALUE.12345.root`

The output file will contain a `TDirectory` named **toys**, which will be empty if no toys are generated (see below for details) and a `TTree` called **limit** with the following branches;

| **Branch name** | **Type** | **Description** |
|------------------------|---------------| ------------------------------------------------------------------------------------------------------------------------|
| **`limit`** | `Double_t` | Main result of combine run with method dependent meaning |
| **`limitErr`** | `Double_t` | Estimated uncertainty on the result |
| **`mh`** | `Double_t` | Value of **MH**, specified with `-m` option |
| **`iToy`** | `Int_t` | Toy number identifier if running with `-t`|
| **`iSeed`** | `Int_t` | Seed specified with `-s`|
| **`t_cpu`** | `Float_t` | Estimated CPU time for algorithm|
| **`t_real`** | `Float_t` | Estimated real time for algorithm|
| **`quantileExpected`** | `Float_t` | Quantile identifier for methods which calculated expected (quantiles) and observed results (eg conversions from $\Delta\ln L$ values) with method dependent meaning. Negative values are reserved for entries which *do not* related to quantiles of a calculation with the default being set to -1 (usually meaning the *observed* result). |

The value of any user defined keyword **$WORD** which is set using `keyword-value` described above will also be included as a  branch with type `string` named **WORD**. The option can be repeated multiple times for multiple keywords.

In some cases, the precise meanings of the branches will depend on the Method being used, which is included in this documentation.

## Toy data generation

By default, each of these methods will be run using the **observed data** as the input. In several cases (as detailed below), it might be useful to run the tool using toy datasets, including Asimov data.

The option `-t` is used to specify to combine to first generate a toy dataset(s) which will be used in replacement of the real data. There are two versions of this,

   * `-t N` with N > 0. Combine will generate N toy datasets from the model and re-run the method once per toy. The seed for the toy generation can be modified with the option `-s` (use `-s -1` for a random seed). The output file will contain one entry in the tree for each of these toys.

   * `-t -1` will produce an Asimov dataset in which statistical fluctuations are suppressed. The procedure to generate this Asimov dataset depends on which type of analysis you are using, see below for details.

The output file will contain the toys (as `RooDataSets` for the observables, including global observables) in the **toys** directory if the option `--saveToys` is provided. If you include this option, the `limit` TTree in the output will have an entry corresponding to the state of the POI used for the generation of the toy, with the value of **`quantileExpected`** set to **-2**. You can add additional branches using the `--trackParameters` option as described in the [common command line options](#common-command-line-options) section above.

!!! warning
    The default values of the nuisance parameters (or any parameter) are used to generate the toy. This means that if, for example, you are using parametric shapes and the parameters inside the workspace are set to arbitrary values, *those* arbitrary values will be used to generate the toy. This behaviour can be modified through the use of the option `--setParameters x=value_x,y=value_y...` which will set the values of the parameters (`x` and `y`) before toy generation. You can also load a snap-shot from a previous fit to set the nuisances to their *post-fit* values (see below).

#### Asimov datasets

If you are using wither `-t -1` or using `AsymptoticLimits`, combine will calculate results based on an Asimov dataset.

   * For counting experiments, the Asimov data will just be set to the total number of expected events (given the values of the nuisance parameters and POIs of the model)

   * For shape analyses with templates, the Asimov dataset will be constructed as a histogram using the same binning which is defined for your analysis.

   * If your model uses parametric shapes (for example when you are using binned data, there are some options as to what Asimov dataset to produce. By *default*, combine will produce the Asimov data as a histogram using the binning which is associated to each observable (ie as set using `RooRealVar::setBins`). If this binning doesn't exist, combine will **guess** a suitable binning - it is therefore best to use `RooRealVar::setBins` to associate a binning to each observable, even if your data is unbinned, if you intend to use Asimov datasets.

You can also ask combine to use a **Pseudo-Asimov** dataset, which is created from many weighted unbinned events.

Setting `--X-rtd TMCSO_AdaptivePseudoAsimov=`$\beta$ with $\beta>0$ will trigger the internal logic of whether to produce a Pseudo-Asimov dataset. This logic is as follows;

   1. For each observable in your dataset, the number of bins, $n_{b}$ is determined either from the value of `RooRealVar::getBins` if it exists or assumed to be 100.

   2. If $N_{b}=\prod_{b}n_{b}>5000$, the number of expected events $N_{ev}$ is determined. *Note* if you are combining multiple channels, $N_{ev}$ refers to the number of expected events in a single channel, the logic is separate for each channel. If  $N_{ev}/N_{b}<0.01$ then a Pseudo-Asimov dataset is created with the number of events equal to $\beta \cdot \mathrm{max}\{100*N_{ev},1000\}$. If $N_{ev}/N_{b}\geq 0.01$ , then a normal Asimov dataset is produced.

   3. If $N_{b}\leq 5000$ then a normal Asimov dataset will be produced

The production of a Pseudo-Asimov dataset can be *forced* by using the option `--X-rtd TMCSO_PseudoAsimov=X` where `X>0` will determine the number of weighted events for the Pseudo-Asimov dataset. You should try different values of `X` since larger values leads to more events in the Pseudo-Asimov dataset resulting in higher precision but in general the fit will be slower. 

You can turn off the internal logic by setting `--X-rtd TMCSO_AdaptivePseudoAsimov=0 --X-rtd TMCSO_PseudoAsimov=0` thereby forcing histograms to be generated.

!!! info
    If you set `--X-rtd TMCSO_PseudoAsimov=X` with `X>0` and also turn on `--X-rtd TMCSO_AdaptivePseudoAsimov=`$\beta$, with $\beta>0$, the internal logic will be used but this time the default will be to generate Pseudo-Asimov datasets, rather than the normal Asimov ones.


#### Nuisance parameter generation

The default method of dealing with systematics is to generate random values (around their nominal values, see above) for the nuisance parameters, according to their prior pdfs centred around their default values, *before* generating the data. The *unconstrained* nuisance parameters (eg `flatParam` or `rateParam`) or those with *flat* priors are **not** randomised before the data generation.

The following are options which define how the toys will be generated,

   * `--toysNoSystematics` the nuisance parameters in each toy are *not* randomised when generating the toy datasets - i.e their nominal values are used to generate the data. Note that for methods which profile (fit) the nuisances, the parameters are still floating when evaluating the likelihood.

   * `--toysFrequentist` the nuisance parameters in each toy are set to their nominal values which are obtained *after fitting first to the data*, with POIs fixed, before generating the data. For evaluating likelihoods, the constraint terms are instead randomised within their Gaussian constraint pdfs around the post-fit nuisance parameter values.

If you are using `toysFrequentist`, be aware that the values set by `--setParameters` will be *ignored* for the toy generation as the *post-fit* values will instead be used (except for any parameter which is also a parameter of interest). You can override this behaviour and choose the nominal values for toy generation for any parameter by adding the option `--bypassFrequentistFit` which will skip the initial fit to data or by loading a snapshot (see below).

!!! warning
    The methods such as `AsymptoticLimits` and `HybridNew --LHCmode LHC-limits`, the  "nominal" nuisance parameter values are taken from fits to the data and are therefore not "blind" to the observed data by default (following the fully frequentist paradigm). See the detailed documentation on these methods for avoiding this and running in a completely "blind" mode.

#### Generate only

It is also possible to generate the toys first and then feed them to the Methods in combine. This can be done using `-M GenerateOnly --saveToys`. The toys can then be read and used with the other methods by specifying `--toysFile=higgsCombineTest.GenerateOnly...` and using the same options for the toy generation.

!!! warning
    Some Methods also use toys within the method itself (eg `AsymptoticLimits` and `HybridNew`). For these, you should **not** specify the toy generation with `-t` or the options above and instead follow the specific instructions.

#### Loading snapshots

Snapshots from workspaces can be loaded and used in order to generate toys using the option `--snapshotName <name of snapshot>`. This will first set the parameters to the values in the snapshot *before* any other parameter options are set and toys are generated.

See the section on [saving post-fit workspaces](/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/#using-best-fit-snapshots) for creating workspaces with post-fit snapshots from `MultiDimFit`.

Here are a few examples of calculations with toys from post-fit workspaces using a workspace with $r, m_{H}$ as parameters of interest

- Throw post-fit toy with b from s+b(floating $r,m_{H}$) fit, s with **r=1.0**, **m=best fit MH**, using nuisance values and constraints re-centered on s+b(floating $r,m_{H}$) fit values (aka frequentist post-fit expected) and compute post-fit expected **r** uncertainty profiling **MH**
    `combine higgsCombinemumhfit.MultiDimFit.mH125.root --snapshotName MultiDimFit -M MultiDimFit --verbose 9 -n randomtest --toysFrequentist --bypassFrequentistFit -t -1 --expectSignal=1 -P r --floatOtherPOIs=1 --algo singles`

- Throw post-fit toy with b from s+b(floating $r,m_{H}$) fit, s with **r=1.0, m=128.0**, using nuisance values and constraints re-centered on s+b(floating $r,m_{H}$) fit values (aka frequentist post-fit expected) and compute post-fit expected significance (with **MH** fixed at 128 implicitly)
    `combine higgsCombinemumhfit.MultiDimFit.mH125.root -m 128 --snapshotName MultiDimFit -M ProfileLikelihood --significance --verbose 9 -n randomtest --toysFrequentist --bypassFrequentistFit --overrideSnapshotMass -t -1 --expectSignal=1 --redefineSignalPOIs r --freezeParameters MH`

- Throw post-fit toy with b from s+b(floating $r,m_{H}$) fit, s with **r=0.0**, using nuisance values and constraints re-centered on s+b(floating $r,m_{H}$) fit values (aka frequentist post-fit expected) and compute post-fit expected and observed asymptotic limit (with **MH** fixed at 128 implicitly)
    `combine higgsCombinemumhfit.MultiDimFit.mH125.root -m 128 --snapshotName MultiDimFit -M AsymptoticLimits --verbose 9 -n randomtest --bypassFrequentistFit --overrideSnapshotMass--redefineSignalPOIs r --freezeParameters MH`

## combineTool for job submission

For longer tasks which cannot be run locally, several methods in combine can be split to run on the *LSF batch* or the *Grid*. The splitting and submission is handled using the `combineTool` (see [this getting started](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#combine-tool) section to get the tool)


### Submission to Condor

The syntax for running on condor with the tool is

```sh
combineTool.py -M ALGO [options] --job-mode condor --sub-opts='CLASSADS' --task-name NAME [--dry-run]
```

with `options` being the usual list of `combine` options. The help option `-h` will give a list of both `combine` and `combineTool` sets of options. This can be used with several different methods from `combine`.

The `--sub-opts` option takes a string with the different ClassAds that you want to set, separated by `\n` as argument (e.g. `'+JobFlavour="espresso"\nRequestCpus=1'`).

The `--dry-run` option will show what will be run without actually doing so / submitting the jobs.

For example, to generate toys (eg for use with limit setting) users running on lxplus at CERN the **condor** mode can be used eg

```sh
combineTool.py -d workspace.root -M HybridNew --LHCmode LHC-limits --clsAcc 0  -T 2000 -s -1 --singlePoint 0.2:2.0:0.05 --saveHybridResult -m 125 --job-mode condor --task-name condor-test --sub-opts='+JobFlavour="tomorrow"'
``` 
The `--singlePoint` option is over-ridden so that this will produce a script for each value of the POI in the range 0.2 to 2.0 in steps of 0.05. You can merge multiple points into a script using `--merge` - e.g adding `--merge 10` to the above command will mean that each job contains *at most* 10 of the values. The scripts are labelled by the `--task-name` option. These will be submitted directly to condor adding any options in `--sub-opts` to the condor submit script. Make sure multiple options are separated by `\n`. The jobs will run and produce output in the **current directory**.

Below is an example for splitting points in a multi-dimensional likelihood scan.

#### Splitting jobs for a multi-dimensional likelihood scan

The option `--split-points` issues the command to split the jobs for `MultiDimFit` when using `--algo grid`. The following example will split the jobs such that there are **10 points** in each of the jobs, which will be submitted to the **8nh** queue.

```sh
combineTool.py datacard.txt -M MultiDimFit --algo grid --points 50 --rMin 0 --rMax 1 --job-mode condor --split-points 10 --sub-opts='+JobFlavour="workday"' --task-name mytask -n mytask
```

Remember, any usual options (such as redefining POIs or freezing parameters) are passed to combine and can be added to the command line for `combineTool`.

!!! info
    The option `-n NAME` should be included to avoid overwriting output files as the jobs will be run inside the directory from which the command is issued.


### Running combine jobs on the Grid

For more CPU-intensive tasks, for example determining limits for complex models using toys, it is generally not feasible to compute all the results interactively. Instead, these jobs can be submitted to the Grid.

In this example we will use the `HybridNew` method of combine to determine an upper limit for a sub-channel of the Run 1 SM $H\rightarrow\tau\tau$ analysis. For full documentation, see the section on [computing limits with toys](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/#computing-limits-with-toys).

With this model it would take too long to find the limit in one go, so instead we create a set of jobs in which each one throws toys and builds up the test statistic distributions for a fixed value of the signal strength. These jobs can then be submitted to a batch system or to the Grid using `crab3`. From the set of output distributions it is possible to extract the expected and observed limits.

For this we will use `combineTool.py`

First we need to build a workspace from the [$H\rightarrow\tau\tau$ datacard](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/81x-root606-integration/data/tutorials/htt/125/htt_tt.txt),

```sh
$ text2workspace.py data/tutorials/htt/125/htt_mt.txt -m 125
$ mv data/tutorials/htt/125/htt_mt.root ./
```

To get an idea of the range of signal strength values we will need to build test-statistic distributions for we will first use the `AsymptoticLimits` method of combine,

```nohighlight
$ combine -M Asymptotic htt_mt.root -m 125
 << Combine >>
[...]
 -- AsymptoticLimits (CLs) --
Observed Limit: r < 1.7384
Expected  2.5%: r < 0.4394
Expected 16.0%: r < 0.5971
Expected 50.0%: r < 0.8555
Expected 84.0%: r < 1.2340
Expected 97.5%: r < 1.7200
```

Based on this, a range of 0.2 to 2.0 should be suitable.

We can use the same command for generating the distribution of test statistics with `combineTool`. The `--singlePoint` option is now enhanced to support expressions that generate a set of calls to combine with different values. The accepted syntax is of the form **MIN:MAX:STEPSIZE**, and multiple comma-separated expressions can be specified.

The script also adds an option `--dry-run` which will not actually call combine but just prints out the commands that would be run, e.g,

```sh
combineTool.py -M HybridNew -d htt_mt.root --LHCmode LHC-limits --singlePoint 0.2:2.0:0.2 -T 2000 -s -1 --saveToys --saveHybridResult -m 125 --dry-run
...
[DRY-RUN]: combine -d htt_mt.root --LHCmode LHC-limits -T 2000 -s -1 --saveToys --saveHybridResult -M HybridNew -m 125 --singlePoint 0.2 -n .Test.POINT.0.2
[DRY-RUN]: combine -d htt_mt.root --LHCmode LHC-limits -T 2000 -s -1 --saveToys --saveHybridResult -M HybridNew -m 125 --singlePoint 0.4 -n .Test.POINT.0.4
[...]
[DRY-RUN]: combine -d htt_mt.root --LHCmode LHC-limits -T 2000 -s -1 --saveToys --saveHybridResult -M HybridNew -m 125 --singlePoint 2.0 -n .Test.POINT.2.0
```

When the `--dry-run` option is removed each command will be run in sequence.

#### Grid submission

Submission to the grid with `crab3` works in a similar way. Before doing so ensure that the `crab3` environment has been sourced, then for compatibility reasons source the CMSSW environment again. We will use the example of generating a grid of test-statistic distributions for limits.

```sh
$ source /cvmfs/cms.cern.ch/crab3/crab.sh; cmsenv
$ combineTool.py -d htt_mt.root -M HybridNew --LHCmode LHC-limits --clsAcc 0 -T 2000 -s -1 --singlePoint 0.2:2.0:0.05 --saveToys --saveHybridResult -m 125 --job-mode crab3 --task-name grid-test --custom-crab custom_crab.py
```

The option `--custom-crab` should point to a python file python containing a function of the form `custom_crab(config)` that will be used to modify the default crab configuration. You can use this to set the output site to your local grid site, or modify other options such as the voRole, or the site blacklist/whitelist.

For example

```python
def custom_crab(config):
  print '>> Customising the crab config'
  config.Site.storageSite = 'T2_CH_CERN'
  config.Site.blacklist = ['SOME_SITE', 'SOME_OTHER_SITE']
```

Again it is possible to use the option `--dry-run` to see what the complete crab config will look like before actually submitted it.

Once submitted the progress can be monitored using the standard `crab` commands. When all jobs are completed copy the output from your sites storage element to the local output folder.

```sh
$ crab getoutput -d crab_grid-test
# Now we have to un-tar the output files
$ cd crab_grid-test/results/
$ for f in *.tar; do tar xf $f; done
$ mv higgsCombine*.root ../../
$ cd ../../
```

These output files should be combined with `hadd`, after which we invoke combine as usual to calculate observed and expected limits from the merged grid as usual.

