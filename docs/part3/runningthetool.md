# How to run the tool

The executable <span style="font-variant:small-caps;">Combine</span> provided by the package is used to invoke the tools via the command line. The statistical analysis method, as well as user settings, are also specified on the command line. To see the full list of available options, you can run:

```sh
combine --help
```

The option `-M` is used to choose the statistical evaluation method. There are several groups of statistical methods:

-   **Asymptotic** likelihood methods:
    -   `AsymptoticLimits`: limits calculated according to the asymptotic formulae in [arxiv:1007.1727](http://arxiv.org/abs/1007.1727).
    -   `Significance`: simple profile likelihood approximation, for calculating significances.
-   **Bayesian** methods:
    -   `BayesianSimple`: performing a classical numerical integration (for simple models only).
    -   `MarkovChainMC`: performing Markov Chain integration, for arbitrarily complex models.
-   **Frequentist** or hybrid bayesian-frequentist methods:
    -   `HybridNew`: compute modified frequentist limits, significance/p-values and confidence intervals according to several possible prescriptions with toys. 
-   **Fitting**
    - `FitDiagnostics`: performs maximum likelihood fits to extract the signal rate, and provides diagnostic tools such as pre- and post-fit figures and correlations
    - `MultiDimFit`: performs maximum likelihood fits and likelihood scans with an arbitrary number of parameters of interest.
-   **Miscellaneous** other modules that do not compute limits or confidence intervals, but use the same framework:
    - `GoodnessOfFit`: perform a goodness of fit test for models including shape information. Several GoF tests are implemented.
    - `ChannelConsistencyCheck`: study the consistency between individual channels in a combination.
    - `GenerateOnly`: generate random or asimov toy datasets for use as input to other methods

The command help is organized into five parts:

-   The *Main options* section indicates how to pass the datacard as input to the tool (`-d datacardName`), how to choose the statistical method (`-M MethodName`), and how to set the verbosity level `-v`
-   Under *Common statistics options*, options common to different statistical methods are given. Examples are `--cl`, to specify the confidence level (default is 0.95), or `-t`, to give the number of toy MC extractions required.
-   The *Common input-output options* section includes, for example, the options to specify the mass hypothesis under study (`-m`) or to include a specific string in the output filename (`--name`). 
-   *Common miscellaneous options*.
-   Further method-specific options are available for each method. By passing the method name via the `-M` option, along with `--help`, the options for that specific method are shown in addition to the common options. 

Not all the available options are discussed in this online documentation; use `--help` to get the documentation of all options.

## Common command-line options

There are a number of useful command-line options that can be used to alter the model (or parameters of the model) at run time. The most commonly used, generic options, are:

-   `-H`: first run a different, faster, algorithm (e.g. the `ProfileLikelihood` described below) to obtain an approximate indication of the limit, which will allow the precise chosen algorithm to converge more quickly. We **strongly recommend** to use this option when using the `MarkovChainMC`, `HybridNew` or `FeldmanCousins` calculators, unless you know in which range your limit lies and you set this range manually (the default is `[0, 20]`)

-   `--rMax`, `--rMin`: manually restrict the range of signal strengths to consider. For Bayesian limits with MCMC, a rule of thumb is that `rMax` should be 3-5 times the limit (a too small value of `rMax` will bias your limit towards low values, since you are restricting the integration range, while a too large value will bias you to higher limits)

-   `--setParameters name=value[,name2=value2,...]` sets the starting values of the parameters, useful e.g. when generating toy MC or when setting the parameters as fixed. This option supports the use of regular expressions by replacing `name` with `rgx{some regular expression}`.

-   `--setParameterRanges name=min,max[:name2=min2,max2:...]` sets the ranges of the parameters (useful e.g. for scans in `MultiDimFit`, or for Bayesian integration). This option supports the use of regular expressions by replacing `name` with `rgx{some regular expression}`.

-   `--redefineSignalPOIs name[,name2,...]` redefines the set of parameters of interest.
    -   If the parameters were constant in the input workspace, they are set to be floating.
    -   Nuisance parameters promoted to parameters of interest are removed from the list of nuisances, and thus they are not randomized in methods that randomize nuisances (e.g. `HybridNew` in non-frequentist mode, or `BayesianToyMC`, or in toy generation with `-t` but without `--toysFreq`). This does not have any impact on algorithms that do not randomize nuisance parameters (e.g. fits, `AsymptoticLimits`, or `HybridNew` in fequentist mode) or on algorithms that treat all parameters in the same way (e.g. `MarkovChainMC`).
    -   Note that constraint terms for the nuisances are **dropped** after promotion to a POI using `--redefineSignalPOI`. To produce a likelihood scan for a nuisance parameter, using `MultiDimFit` with **`--algo grid`**, you should instead use the `--parameters (-P)` option, which will not cause the loss of the constraint term when scanning.
    -   Parameters of interest of the input workspace that are not selected by this command become unconstrained nuisance parameters, but they are not added to the list of nuisances so they will not be randomized (see above).


-   `--freezeParameters name1[,name2,...]` Will freeze the parameters with the given names to their set values. This option supports the use of regular expression by replacing `name` with `rgx{some regular expression}` for matching to *constrained nuisance parameters* or `var{some regular expression}` for matching to *any* parameter. For example `--freezeParameters rgx{CMS_scale_j.*}` will freeze all constrained nuisance parameters with the prefix `CMS_scale_j`, while `--freezeParameters var{.*rate_scale}` will freeze any parameter (constrained nuisance parameter or otherwise) with the suffix `rate_scale`.
    - Use the option `--freezeParameters allConstrainedNuisances` to freeze all nuisance parameters that have a constraint term (i.e not `flatParams` or `rateParams` or other freely floating parameters).
    - Similarly, the option `--floatParameters name1[,name2,...]` sets the parameter(s) floating and also accepts regular expressions.
    - Groups of nuisance parameters (constrained or otherwise), as defined in the datacard, can be frozen using `--freezeNuisanceGroups`. You can also freeze all nuisances that are *not* contained in a particular group using a **^** before the group name (`--freezeNuisanceGroups=^group_name` will freeze everything except nuisance parameters in the group "group_name".)
    - All *constrained* nuisance parameters (not `flatParam` or `rateParam`) can be set floating using `--floatAllNuisances`.

!!! warning
    Note that the floating/freezing options have a priority ordering from lowest to highest as `floatParameters < freezeParameters < freezeNuisanceGroups < floatAllNuisances`. Options with higher priority will take precedence over those with lower priority.

-   `--trackParameters name1[,name2,...]` will add a branch to the output tree for each of the named parameters. This option supports the use of regular expressions by replacing `name` with `rgx{some regular expression}`

    - The name of the branch will be **trackedParam_*name***.
    - The exact behaviour depends on the method used. For example, when using `MultiDimFit` with `--algo scan`, the value of the parameter at each point in the scan will be saved, while for `FitDiagnostics`, only the value at the end of the fit will be saved.

-   `--trackErrors name1[,name2,...]` will add a branch to the output tree for the error of each of the named parameters. This option supports the use of regular expressions by replacing `name` with `rgx{some regular expression}`

    - The name of the branch will be **trackedError_*name***.
    - The behaviour, in terms of which values are saved, is the same as `--trackParameters` above.

By default, the data set used by <span style="font-variant:small-caps;">Combine</span> will be the one listed in the datacard. You can tell <span style="font-variant:small-caps;">Combine</span> to use a different data set (for example a toy data set that you generated) by using the option `--dataset`. The argument should be `rootfile.root:workspace:location` or `rootfile.root:location`. In order to use this option, you must first convert your datacard to a binary workspace and use this binary workspace as the input to <span style="font-variant:small-caps;">Combine</span>. 

### Generic Minimizer Options

<span style="font-variant:small-caps;">Combine</span> uses its own minimizer class, which is used to steer Minuit (via RooMinimizer), named the `CascadeMinimizer`. This allows for sequential minimization, which can help in case a particular setting or algorithm fails. The `CascadeMinimizer` also knows about extra features of <span style="font-variant:small-caps;">Combine</span> such as *discrete* nuisance parameters.

All of the fits that are performed in <span style="font-variant:small-caps;">Combine</span>'s methods use this minimizer. This means that the fits can be tuned using these common options,

* `--cminPoiOnlyFit`: First, perform a fit floating *only* the parameters of interest. This can be useful to find, roughly, where the global minimum is.
* `--cminPreScan`: Do a scan before the first minimization.
* `--cminPreFit arg` If set to a value N > 0, the minimizer will perform a pre-fit with strategy (N-1), with the nuisance parameters frozen.
     * `--cminApproxPreFitTolerance arg`: If non-zero, first do a pre-fit with this tolerance (or 10 times the final tolerance, whichever is largest)
     * `--cminApproxPreFitStrategy arg`:   Strategy to use in the pre-fit. The default is strategy 0.
* `--cminDefaultMinimizerType arg`: Set the default minimizer type. By default this is set to Minuit2.
* `--cminDefaultMinimizerAlgo arg`: Set the default minimizer algorithm. The default algorithm is Migrad.
* `--cminDefaultMinimizerTolerance arg`: Set the default minimizer tolerance, the default is 0.1.
* `--cminDefaultMinimizerStrategy arg`: Set the default minimizer strategy between 0 (speed), 1 (balance - *default*), 2 (robustness). The [Minuit documentation](http://www.fresco.org.uk/minuit/cern/node6.html) for this is pretty sparse but in general, 0 means evaluate the function less often, while 2 will waste function calls to get precise answers. An important note is that the `Hesse` algorithm (for error and correlation estimation) will be run *only* if the strategy is 1 or 2.
* `--cminFallbackAlgo arg`: Provides a list of fallback algorithms, to be used in case the default minimizer fails. You can provide multiple options using the syntax `Type[,algo],strategy[:tolerance]`: eg `--cminFallbackAlgo Minuit2,Simplex,0:0.1` will fall back to the simplex algorithm of Minuit2 with strategy 0 and a tolerance 0.1, while `--cminFallbackAlgo Minuit2,1` will use the default algorithm (Migrad) of Minuit2 with strategy 1.
* `--cminSetZeroPoint (0/1)`: Set the reference of the NLL to 0 when minimizing, this can help faster convergence to the minimum if the NLL itself is large. The default is true (1), set to 0 to turn off.

The allowed combinations of minimizer types and minimizer algorithms are as follows:

| **Minimizer type** | **Minimizer algorithm** |
|--------------------|--------------------|
|`Minuit`	      | `Migrad`, `Simplex`, `Combined`, `Scan` |
|`Minuit2` 	      | `Migrad`, `Simplex`, `Combined`, `Scan` |
|`GSLMultiMin`       | `ConjugateFR`, `ConjugatePR`, `BFGS`, `BFGS2`, `SteepestDescent`|

You can find details about these in the Minuit2 documentation [here](https://root.cern.ch/root/htmldoc/guides/minuit2/Minuit2.html).

More of these options can be found in the **Cascade Minimizer options** section when running `--help`.


### Output from combine

Most methods will print the results of the computation to the screen. However, in addition, <span style="font-variant:small-caps;">Combine</span> will also produce a root file containing a tree called **limit** with these results. The name of this file will be of the format,

    higgsCombineTest.MethodName.mH$MASS.[word$WORD].root

where **$WORD** is any user defined keyword from the datacard which has been set to a particular value.

A few command-line options can be used to control this output:

-   The option `-n` allows you to specify part of the name of the root file. e.g. if you pass `-n HWW` the root file will be called `higgsCombineHWW....` instead of `higgsCombineTest`
-   The option `-m` allows you to specify the (Higgs boson) mass hypothesis, which gets written in the filename and in the output tree. This simplifies the bookeeping, as it becomes possible to merge multiple trees corresponding to different (Higgs boson) masses using `hadd`. Quantities can then be plotted as a function of the mass. The default value is m=120.
-   The option `-s` can be used to specify the seed (eg `-s 12345`) used in toy generation. If this option is given, the name of the file will be extended by this seed, eg `higgsCombineTest.AsymptoticLimits.mH120.12345.root`
-   The option `--keyword-value` allows you to specify the value of a keyword in the datacard such that **$WORD** (in the datacard) will be given the value of **VALUE** in the command `--keyword-value WORD=VALUE`, eg  `higgsCombineTest.AsymptoticLimits.mH120.WORDVALUE.12345.root`

The output file will contain a `TDirectory` named **toys**, which will be empty if no toys are generated (see below for details) and a `TTree` called **limit** with the following branches;

| **Branch name** | **Type** | **Description** |
|------------------------|---------------| ------------------------------------------------------------------------------------------------------------------------|
| **`limit`** | `Double_t` | Main result of combine run, with method-dependent meaning |
| **`limitErr`** | `Double_t` | Estimated uncertainty on the result |
| **`mh`** | `Double_t` | Value of **MH**, specified with `-m` option |
| **`iToy`** | `Int_t` | Toy number identifier if running with `-t`|
| **`iSeed`** | `Int_t` | Seed specified with `-s`|
| **`t_cpu`** | `Float_t` | Estimated CPU time for algorithm|
| **`t_real`** | `Float_t` | Estimated real time for algorithm|
| **`quantileExpected`** | `Float_t` | Quantile identifier for methods that calculated expected (quantiles) and observed results (eg conversions from $\Delta\ln L$ values), with method-dependent meaning. Negative values are reserved for entries that *do not* relate to quantiles of a calculation, with the default being set to -1 (usually meaning the *observed* result). |

The value of any user-defined keyword **$WORD** that is set using `keyword-value` described above will also be included as a branch with type `string` named **WORD**. The option can be repeated multiple times for multiple keywords.

In some cases, the precise meanings of the branches will depend on the method being used. In this case, it will be specified in this documentation.

## Toy data generation

By default, each of the methods described so far will be run using the **observed data** as the input. In several cases (as detailed below), it is useful to run the tool using toy datasets, including Asimov data sets.

The option `-t` is used to tell <span style="font-variant:small-caps;">Combine</span> to first generate one or more toy data sets, which will be used instead of the observed data. There are two versions,

   * `-t N` with N > 0. <span style="font-variant:small-caps;">Combine</span> will generate N toy datasets from the model and re-run the method once per toy. The seed for the toy generation can be modified with the option `-s` (use `-s -1` for a random seed). The output file will contain one entry in the tree for each of these toys.

   * `-t -1` will produce an Asimov data set, in which statistical fluctuations are suppressed. The procedure for generating this Asimov data set depends on the type of analysis you are using. More details are given below. 

!!! warning
    The default values of the nuisance parameters (or any parameter) are used to generate the toy. This means that if, for example, you are using parametric shapes and the parameters inside the workspace are set to arbitrary values, *those* arbitrary values will be used to generate the toy. This behaviour can be modified through the use of the option `--setParameters x=value_x,y=value_y...`, which will set the values of the parameters (`x` and `y`) before toy generation. You can also load a snapshot from a previous fit to set the nuisance parameters to their *post-fit* values (see below).

The output file will contain the toys (as `RooDataSets` for the observables, including global observables) in the **toys** directory if the option `--saveToys` is provided. If you include this option, the `limit` TTree in the output will have an entry corresponding to the state of the POI used for the generation of the toy, with the value of **`quantileExpected`** set to **-2**. 

The branches that are created by methods like `MultiDimFit` *will not* show the values used to generate the toy. If you also want the TTree to show the values of the POIs used to generate the toy, you should add additional branches using the `--trackParameters` option as described in the [common command-line options](#common-command-line-options) section above. These branches will behave as expected when adding the option `--saveToys`. 

!!! warning
    For statistical methods that make use of toys (including `HybridNew`, `MarkovChainMC` and running with `-t N`), the results of repeated <span style="font-variant:small-caps;">Combine</span> commands will not be identical when using the datacard as the input. This is due to a feature in the tool that allows one to run concurrent commands that do not interfere with one another. In order to produce reproducible results with toy-based methods, you should first convert the datacard to a binary workspace using `text2workspace.py` and then use the resulting file as input to the <span style="font-variant:small-caps;">Combine</span> commands
    


### Asimov datasets

If you are using either `-t -1` or  `AsymptoticLimits`, <span style="font-variant:small-caps;">Combine</span> will calculate results based on an Asimov data set.

   * For counting experiments, the Asimov data set will just be the total number of expected events (given the values of the nuisance parameters and POIs of the model)

   * For shape analyses with templates, the Asimov data set will be constructed as a histogram using the same binning that is defined for your analysis.

   * If your model uses parametric shapes, there are some options as to what Asimov data set to produce. By *default*, <span style="font-variant:small-caps;">Combine</span> will produce the Asimov data set as a histogram using the binning that is associated with each observable (ie as set using `RooRealVar::setBins`). If this binning does not exist, <span style="font-variant:small-caps;">Combine</span> will **guess** a suitable binning - it is therefore best to use `RooRealVar::setBins` to associate a binning with each observable, even if your data is unbinned, if you intend to use Asimov data sets.

You can also ask <span style="font-variant:small-caps;">Combine</span> to use a **Pseudo-Asimov** dataset, which is created from many weighted unbinned events.

Setting `--X-rtd TMCSO_AdaptivePseudoAsimov=`$\beta$ with $\beta>0$ will trigger the internal logic of whether to produce a Pseudo-Asimov dataset. This logic is as follows;

   1. For each observable in your dataset, the number of bins, $n_{b}$ is determined either from the value of `RooRealVar::getBins`, if it exists, or assumed to be 100.

   2. If $N_{b}=\prod_{b}n_{b}>5000$, the number of expected events $N_{ev}$ is determined. *Note* if you are combining multiple channels, $N_{ev}$ refers to the number of expected events in a single channel. The logic is separate for each channel. If  $N_{ev}/N_{b}<0.01$ then a Pseudo-Asimov data set is created with the number of events equal to $\beta \cdot \mathrm{max}\{100*N_{ev},1000\}$. If $N_{ev}/N_{b}\geq 0.01$ , then a normal Asimov data set is produced.

   3. If $N_{b}\leq 5000$ then a normal Asimov data set will be produced

The production of a Pseudo-Asimov data set can be *forced* by using the option `--X-rtd TMCSO_PseudoAsimov=X` where `X>0` will determine the number of weighted events for the Pseudo-Asimov data set. You should try different values of `X`, since larger values lead to more events in the Pseudo-Asimov data set, resulting in higher precision. However, in general, the fit will be slower. 

You can turn off the internal logic by setting `--X-rtd TMCSO_AdaptivePseudoAsimov=0 --X-rtd TMCSO_PseudoAsimov=0`, thereby forcing histograms to be generated.

!!! info
    If you set `--X-rtd TMCSO_PseudoAsimov=X` with `X>0` and also turn on `--X-rtd TMCSO_AdaptivePseudoAsimov=`$\beta$, with $\beta>0$, the internal logic will be used, but this time the default will be to generate Pseudo-Asimov data sets, rather than the standard Asimov ones.


### Nuisance parameter generation

The default method of handling systematics is to generate random values (around their nominal values, see above) for the nuisance parameters, according to their prior PDFs centred around their default values, *before* generating the data. The *unconstrained* nuisance parameters (eg `flatParam` or `rateParam`), or those with *flat* priors are **not** randomized before the data generation. If you wish to also randomize these parameters, you **must** declare them as `flatParam` in your datacard and, when running text2workspace, you must add the option `--X-assign-flatParam-prior` to the command line.

The following options define how the toys will be generated,

   * `--toysNoSystematics` the nuisance parameters in each toy are *not* randomized when generating the toy data sets - i.e their nominal values are used to generate the data. Note that for methods which profile (fit) the nuisances, the parameters are still floating when evaluating the likelihood.

   * `--toysFrequentist` the nuisance parameters in each toy are set to their nominal values which are obtained *after first fitting to the observed data*, with the POIs fixed, before generating the toy data sets. For evaluating likelihoods, the constraint terms are instead randomized within their PDFs around the post-fit nuisance parameter values.

If you are using `toysFrequentist`, be aware that the values set by `--setParameters` will be *ignored* for the toy generation as the *post-fit* values will instead be used (except for any parameter that is also a parameter of interest). You can override this behaviour and choose the nominal values for toy generation for any parameter by adding the option `--bypassFrequentistFit`, which will skip the initial fit to data, or by loading a snapshot (see below).

!!! warning
    For methods such as `AsymptoticLimits` and `HybridNew --LHCmode LHC-limits`, the  "nominal" nuisance parameter values are taken from fits to the data and are, therefore, not "blind" to the observed data by default (following the fully frequentist paradigm). See the detailed documentation on these methods for how to run in fully "blinded" mode.

### Generate only

It is also possible to generate the toys first, and then feed them to the methods in <span style="font-variant:small-caps;">Combine</span>. This can be done using `-M GenerateOnly --saveToys`. The toys can then be read and used with the other methods by specifying `--toysFile=higgsCombineTest.GenerateOnly...` and using the same options for the toy generation. 

You can specify to run on a single toy, in place of the observed data, by including the option `-D file.root:toys/toy_i`. For example adding `-D higgsCombineTest.GenerateOnly.mH120.123456.root:toys/toy_10` will run on the  data set `toy_10` (the 10th toy) that was generated and saved in the file `higgsCombineTest.GenerateOnly.mH120.123456.root`. 

!!! warning
    Some methods also use toys within the method itself (eg `AsymptoticLimits` and `HybridNew`). For these, you should **not** specify the toy generation with `-t` or the options above. Instead, you should follow the method-specific instructions.

### Loading snapshots

Snapshots from workspaces can be loaded and used in order to generate toys using the option `--snapshotName <name of snapshot>`. This will first set the parameters to the values in the snapshot, *before* any other parameter options are set and toys are generated.

See the section on [saving post-fit workspaces](/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/#using-best-fit-snapshots) for creating workspaces with post-fit snapshots from `MultiDimFit`.

Here are a few examples of calculations with toys from post-fit workspaces using a workspace with $r, m_{H}$ as parameters of interest

- Throw post-fit toy with b from s+b(floating $r,m_{H}$) fit, s with **r=1.0**, **m=best fit MH**, using nuisance parameter values and constraints re-centered on s+b(floating $r,m_{H}$) fit values (aka frequentist post-fit expected) and compute post-fit expected **r** uncertainty profiling **MH**
    `combine higgsCombinemumhfit.MultiDimFit.mH125.root --snapshotName MultiDimFit -M MultiDimFit --verbose 9 -n randomtest --toysFrequentist --bypassFrequentistFit -t -1 --expectSignal=1 -P r --floatOtherPOIs=1 --algo singles`

- Throw post-fit toy with b from s+b(floating $r,m_{H}$) fit, s with **r=1.0, m=128.0**, using nuisance parameter values and constraints re-centered on s+b(floating $r,m_{H}$) fit values (aka frequentist post-fit expected) and compute post-fit expected significance (with **MH** fixed at 128 implicitly)
    `combine higgsCombinemumhfit.MultiDimFit.mH125.root -m 128 --snapshotName MultiDimFit -M ProfileLikelihood --significance --verbose 9 -n randomtest --toysFrequentist --bypassFrequentistFit --overrideSnapshotMass -t -1 --expectSignal=1 --redefineSignalPOIs r --freezeParameters MH`

- Throw post-fit toy with b from s+b(floating $r,m_{H}$) fit, s with **r=0.0**, using nuisance parameter values and constraints re-centered on s+b(floating $r,m_{H}$) fit values (aka frequentist post-fit expected) and compute post-fit expected and observed asymptotic limit (with **MH** fixed at 128 implicitly)
    `combine higgsCombinemumhfit.MultiDimFit.mH125.root -m 128 --snapshotName MultiDimFit -M AsymptoticLimits --verbose 9 -n randomtest --bypassFrequentistFit --overrideSnapshotMass--redefineSignalPOIs r --freezeParameters MH`

## combineTool for job submission

For longer tasks that cannot be run locally, several methods in <span style="font-variant:small-caps;">Combine</span> can be split to run on a *batch* system or on the *Grid*. The splitting and submission is handled using the `combineTool.py` script.


### Submission to Condor

The syntax for running on condor with the tool is

```sh
combineTool.py -M ALGO [options] --job-mode condor --sub-opts='CLASSADS' --task-name NAME [--dry-run]
```

with `options` being the usual list of <span style="font-variant:small-caps;">Combine</span> options. The help option `-h` will give a list of both <span style="font-variant:small-caps;">Combine</span> and `combineTool` options. It is possible to use this tool with several different methods from <span style="font-variant:small-caps;">Combine</span>.

The `--sub-opts` option takes a string with the different ClassAds that you want to set, separated by `\n` as argument (e.g. `'+JobFlavour="espresso"\nRequestCpus=1'`).

The `--dry-run` option will show what will be run without actually doing so / submitting the jobs.

For example, to generate toys (eg for use with limit setting) users running on lxplus at CERN can use the **condor** mode:

```sh
combineTool.py -d workspace.root -M HybridNew --LHCmode LHC-limits --clsAcc 0  -T 2000 -s -1 --singlePoint 0.2:2.0:0.05 --saveHybridResult -m 125 --job-mode condor --task-name condor-test --sub-opts='+JobFlavour="tomorrow"'
``` 
The `--singlePoint` option is over-ridden, so that this will produce a script for each value of the POI in the range 0.2 to 2.0 in steps of 0.05. You can merge multiple points into a script using `--merge` - e.g adding `--merge 10` to the above command will mean that each job contains *at most* 10 of the values. The scripts are labelled by the `--task-name` option. They will be submitted directly to condor, adding any options in `--sub-opts` to the condor submit script. Make sure multiple options are separated by `\n`. The jobs will run and produce output in the **current directory**.

Below is an example for splitting points in a multi-dimensional likelihood scan.

#### Splitting jobs for a multi-dimensional likelihood scan

The option `--split-points` issues the command to split the jobs for `MultiDimFit` when using `--algo grid`. The following example will split the jobs such that there are **10 points** in each of the jobs, which will be submitted to the **workday** queue.

```sh
combineTool.py datacard.txt -M MultiDimFit --algo grid --points 50 --rMin 0 --rMax 1 --job-mode condor --split-points 10 --sub-opts='+JobFlavour="workday"' --task-name mytask -n mytask
```

Remember, any usual options (such as redefining POIs or freezing parameters) are passed to <span style="font-variant:small-caps;">Combine</span> and can be added to the command line for `combineTool`.

!!! info
    The option `-n NAME` should be included to avoid overwriting output files, as the jobs will be run inside the directory from which the command is issued.


### Grid submission with combineTool

For more CPU-intensive tasks, for example determining limits for complex models using toys, it is generally not feasible to compute all the results interactively. Instead, these jobs can be submitted to the Grid.

In this example we will use the `HybridNew` method of <span style="font-variant:small-caps;">Combine</span> to determine an upper limit for a sub-channel of the Run 1 SM $H\rightarrow\tau\tau$ analysis. For full documentation, see the section on [computing limits with toys](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/#computing-limits-with-toys).

With this model it would take too long to find the limit in one go, so instead we create a set of jobs in which each one throws toys and builds up the test statistic distributions for a fixed value of the signal strength. These jobs can then be submitted to a batch system or to the Grid using `crab3`. From the set of output distributions it is possible to extract the expected and observed limits.

For this we will use `combineTool.py`

First we need to build a workspace from the [$H\rightarrow\tau\tau$ datacard](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/data/tutorials/htt/125/htt_tt.txt),

```sh
$ text2workspace.py data/tutorials/htt/125/htt_mt.txt -m 125
$ mv data/tutorials/htt/125/htt_mt.root ./
```

To get an idea of the range of signal strength values we will need to build test-statistic distributions for, we will first use the `AsymptoticLimits` method of <span style="font-variant:small-caps;">Combine</span>,

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

We can use the same command for generating the distribution of test statistics with `combineTool`. The `--singlePoint` option is now enhanced to support expressions that generate a set of calls to <span style="font-variant:small-caps;">Combine</span> with different values. The accepted syntax is of the form **MIN:MAX:STEPSIZE**, and multiple comma-separated expressions can be specified.

The script also adds an option `--dry-run`, which will not actually call com<span style="font-variant:small-caps;">Combine</span>bine but just prints out the commands that would be run, e.g,

```sh
combineTool.py -M HybridNew -d htt_mt.root --LHCmode LHC-limits --singlePoint 0.2:2.0:0.2 -T 2000 -s -1 --saveToys --saveHybridResult -m 125 --dry-run
...
[DRY-RUN]: combine -d htt_mt.root --LHCmode LHC-limits -T 2000 -s -1 --saveToys --saveHybridResult -M HybridNew -m 125 --singlePoint 0.2 -n .Test.POINT.0.2
[DRY-RUN]: combine -d htt_mt.root --LHCmode LHC-limits -T 2000 -s -1 --saveToys --saveHybridResult -M HybridNew -m 125 --singlePoint 0.4 -n .Test.POINT.0.4
[...]
[DRY-RUN]: combine -d htt_mt.root --LHCmode LHC-limits -T 2000 -s -1 --saveToys --saveHybridResult -M HybridNew -m 125 --singlePoint 2.0 -n .Test.POINT.2.0
```

When the `--dry-run` option is removed each command will be run in sequence.

### Grid submission with crab3

Submission to the grid with `crab3` works in a similar way. Before doing so, ensure that the `crab3` environment has been sourced in addition to the CMSSW environment. We will use the example of generating a grid of test-statistic distributions for limits.

```sh
$ cmsenv; source /cvmfs/cms.cern.ch/crab3/crab.sh
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

Again it is possible to use the option `--dry-run` to see what the complete crab config will look like before actually submitting it.

Once submitted, the progress can be monitored using the standard `crab` commands. When all jobs are completed, copy the output from your site's storage element to the local output folder.

```sh
$ crab getoutput -d crab_grid-test
# Now we have to un-tar the output files
$ cd crab_grid-test/results/
$ for f in *.tar; do tar xf $f; done
$ mv higgsCombine*.root ../../
$ cd ../../
```

These output files should be combined with `hadd`, after which we invoke <span style="font-variant:small-caps;">Combine</span> as usual to calculate observed and expected limits from the merged grid, as usual.

