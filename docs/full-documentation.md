Full documentation
------------------

### Common options

#### Choice of the prior for bayesian methods

For methods based on bayesian inference (BayesianSimple, MarkovChainMC) you can specify a prior for the signal strength. This is controlled through option **`--prior`**, and the values that are currently allowed are:

-   **uniform**: flat prior
-   **1/sqrt(r)**: prior proportional to *1/sqrt(r)*, where *r* is the signal strength; for a counting experiment with a single channel and no systematics, this is Jeffrey's prior.
    Note that you might have to put the *1/sqrt(r)* within quotes because for some shells the braces have a special meaning.

If no prior is specified, a flat prior is assumed.

### Algorithm-specific options

\#ProfiLe

#### ProfileLikelihood algorithm

The ProfileLikelihood algorithm has only two options besides the common ones: the choice of the minimizer algorith (Minuit2 or Minuit), and the choice of the tolerance.

If you see that the limit fails with one of the minimizer, try the other. Sometimes this problem happens due to numerical instabilities even if the model itself looks perfectly well behaved. If neither of the two succeeds, another possible way of circumventing the instability is sometimes to change the range of `r` by a small amount, or to change slightly one number in the model (e.g. in one simple counting experiment a test was failing when assuming ΔB/B = 0.3, but was succeeding for ΔB/B = 0.301 and ΔB/B = 0.299, giving almost the same result)

In case you experience numerical instabilities, e.g. failures or bogus results, you could be able to work around the problem by performing the minimization multiple times and requiring the results to be consistent. This can be done using the options below.

-   **`maxTries`**: how many times the program should attempt to find a minimum at most (plausible range: 20-100)
-   **`tries`**: how many times the program should succeed in computing a limit (plausible range: 5-10)
-   **`maxOutlierFraction`**: maximum fraction of bad results within the tries ones (plausible: 0.15-0.30; default 0.25)
-   **`maxOutliers`**: maximum number of bogus results before giving up for this point (plausible values: 2-5; default 3)
-   **`preFit`**: don't try to profile the likelihood if Minuit already gave errors in finding the minimum. (suggested)

#### BayesianSimple algorithm

The BayesianSimple algorithm doesn't have parameters besides the ones specified above under the "Common options" section.

\#MarkoV

#### MarkovChainMC algorithm

This algorithm has many configurable parameters, and you're encouraged to experiment with those because the default configuration might not be the best for you (or might not even work for you at all)

\*Iterations, burn-in, tries\*
Three parameters control how the MCMC integration is performed:

-   the number of **tries** (option **`--tries`**): the algorithm will run multiple times with different ransom seeds and report as result the truncated mean and rms of the different results. The default value is 10, which should be ok for a quick computation, but for something more accurate you might want to increase this number even up to ~200.
-   the number of **iterations** (option **`-i`**) determines how many points are proposed to fill a single Markov Chain. The default value is 10k, and a plausible range is between 5k (for quick checks) and 20-30k for lengthy calculations. Usually beyond 30k you get a better tradeoff in time vs accuracy by increasing the number of chains (option **`tries`**)
-   the number of **burn-in steps** (option **`-b`**) is the number of points that are removed from the beginning of the chain before using it to compute the limit. IThe default is 200. If your chain is very long, you might want to try increase this a bit (e.g. to some hundreds). Instead going below 50 is probably dangerous.

**Proposal**
The option **`--proposal`** controls the way new points are proposed to fill in the MC chain.

-   **`uniform`**: pick points at random. This works well if you have very few nuisance parameters (or none at all), but normally fails if you have many.
-   **`gaus`**: Use a product of independent gaussians one for each nuisance parameter; the sigma of the gaussian for each variable is 1/5 of the range of the variable (this can be controlled using the parameter **`propHelperWidthRangeDivisor`**). This proposal appears to work well for a reasonable number of nuisances (up to ~15), provided that the range of the nuisance parameters is reasonable, like ±5σ. It does **not** work without systematics.
-   `=ortho` (**default**): (previously known as `test`) This proposalis similar to the multi-gaussian proposal but at every step only a single coordinate of the point is varied, so that the acceptance of the chain is high even for a large number of nuisances (i.e. more than 20).
-   **`fit`**: Run a fit and use the uncertainty matrix from HESSE to construct a proposal (or the one from MINOS if the option `--runMinos` is specified). This sometimes work fine, but sometimes gives biased results, so we don't recommend it in general.

If you believe there's something going wrong, e.g. if your chain remains stuck after accepting only a few events, the option **`debugProposal`** can be used to have a printout of the first *N* proposed points to see what's going on (e.g. if you have some region of the phase space with probability zero, the `gaus` and `fit` proposal can get stuck there forever)

\#HybridNew

#### HybridNew algorithm

**Type of limit**
By default, HybridNew computes hybrid bayesian-frequentist limits. If one specifies the command line option **`--freqentist`** then it will instead compute the fully frequentist limits.

**Rule** (option **`--rule`**)
The rule defines how the distribution of the test statistics is used to compute a limit. When using the *CL&lt;sub&gt;s+b&lt;/sub&gt;* rule the limit to the value of the signal strength for which 95% of the pseudo-experiments give a result more signal-like than the current one, &lt;nobr&gt;*P(x &lt; x&lt;sub&gt;obs&lt;/sub&gt;|r\*s+b) = 1 - CL*&lt;/nobr&gt;. For the more conservative *CL&lt;sub&gt;s&lt;/sub&gt;* rule, the limit is set by the condition &lt;nobr&gt;*P(x &lt; x&lt;sub&gt;obs&lt;/sub&gt;|r\*s+b) /P(x &lt; x&lt;sub&gt;obs&lt;/sub&gt;|b) = 1 - CL*&lt;/nobr&gt; .
**The default rule is *CL&lt;sub&gt;s&lt;/sub&gt;***.

**Test statistics** (option **`--testStat`**)
The **test statistics** is the measure of how signal-like or background-like is an observation. The following test statistics are provided:

-   **Simple Likelihood Ratio** (option value **`LEP`** or **`SLR`**): The LEP-like ratio of likelihoods ( *L(x|r\*s+b,θ) / L(x|b, θ)* ) where numerator and denominator are evaluated for the reference values of the nuisance parameters θ.
-   **Ratio of Profiled Likelihoods** (option value **`TEV`** or **`ROPL`**): The Tevatron-like ratio of profiled likelihoods, in which before computing each of the likelihoods is maximised as function of the nuisance parameters ( *max&lt;sub&gt;θ&lt;/sub&gt; L(x|r\*s+b,θ) / max&lt;sub&gt;θ&lt;/sub&gt; L(x|b, θ)* ).
-   **Profile Likelihood modified for upper limits** (option value **`LHC`** or **`MPL`**): The LHC-like (or Atlas-like) profiled likelihood in which the maximization of the likelihood is done also in the signal strength (*max&lt;sub&gt;θ&lt;/sub&gt; L(x|r\*s+b,θ) / max&lt;sub&gt;r', θ&lt;/sub&gt; L(x|r'\*s+b,θ)* ), with the constraints *0 ≤ r' ≤ r* where the upper bound is applied to force the method to always give an upper limit and not a two-sided interval.
-   **Profile Likelihood (not modified for upper limits)** (option value **`PL`**): The traditional profiled likelihood in which the maximization of the likelihood is done also in the signal strength (*max&lt;sub&gt;θ&lt;/sub&gt;L(x|r\*s+b,θ) / max&lt;sub&gt;r', θ&lt;/sub&gt;L(x|x|r'\*s+b,θ)* ), with just the physical constraints *0 ≤ r'* This test statistics can give two-sided limits, as it starts decreasing when the number of observed events is above the expectations from the signal+background hypothesis.

The default value when computing hybrid bayesian-frequentist limits is **`LEP`**. The default value when computing frequentist limits is **`LHC`**.

**Accuracy**
The search for the limit is performed using an adaptive algorithm, terminating when the estimate of the limit value is below some limit or when the precision cannot be futher improved with the specified options. The options controlling this behaviour are:

-   **`rAbsAcc`**, **`rRelAcc`**: define the accuracy on the limit at which the search stops. The default values are 0.1 and 0.05 respectively, meaning that the search is stopped when Δr &lt; 0.1 or Δr/r &lt; 0.05.
-   **`clsAcc`**: this determines the absolute accuracy up to which the CLs values are computed when searching for the limit. The default is 0.5%. Raising the accuracy above this value will increase significantly the time to run the algorithm, as you need N&lt;sup&gt;2&lt;/sup&gt; more toys to improve the accuracy by a factor &lt;sub&gt;N&lt;/sub&gt;; you can consider enlarging this value if you're computing limits with a larger CL (e.g. 90% or 68%)
    Note that if you're using the *CL&lt;sub&gt;s+b&lt;/sub&gt;* rule then this parameter will control the uncertainty on *CL&lt;sub&gt;s+b&lt;/sub&gt;*, not on *CL&lt;sub&gt;s&lt;/sub&gt;* .
-   **`T`** or **`toysH`**: controls the minimum number of toys that are generated for each point. The default value of 500 should be ok when computing the limit with 90-95% CL. You can decrease this number if you're computing limits at 68% CL, or increase it if you're using 99% CL.

**Computing significances**
When computing the significances, there is no adaptive generation. You can control the number of toys with option **`T`** or `toysH=` and the option **`iterations`** (shortcut **`-i`**, default 1): the default of (1 iteration)×(500 toys) is not enough to probe a significances above ~2. We suggest that you uncrease the number of iterations instead of the number of toys, since the increase in time is linear with the iterations but non-linear with the toys.

In order to compute the significance in multiple jobs, proceed as follows:

-   Run N different jobs with the same inputs but different random seed (option **`-s`**), specifying the additional option **`--saveHybridResult`**.
-   Use **`hadd`** to merge the output root files in a single one. The program will complain with messages like `Cannot merge object type, name: HybridCalculator _result` which can be safely ignored.
-   Compute the significance from the merged file running again but with options `--readHybridResults=` and **`--toysFile=<merged.root>`**

Caveat: there is no check within the program that you're using consistent inputs for the different jobs.

**Simple hypotesis testing**:
Sometimes you just want to compute the p-values of the background-only hypotesis and of the signal plus background hypotesis for a fixed value of the signal strength.
This can be done specifying the option **`singlePoint <value>`** which will set the signal strength to that value and run the hypothesis test. It will generate toys until the required accuracy is met (see above for parameter **`clsAcc`**). You can turn off adaptive generation setting **`clsAcc`** to zero, and then it will generate the toys once (or multiple times if you set option **`iterations`** to a value larger than 1).
Just like for significance, you can run multiple times with different seeds and options **`--saveHybridResult`**, combine the output files with **`hadd`** and then compute the final result with **`--readHybridResult --toysFile=merged.root`**

**Performance issues**
The current hybrid code requires a lot of cpu resources. You can speed up the processing by using multiple cores (option **`fork`**, default value is 1). Note that even with **`fork`** set to 1, toy generation is done in a separate thread to avoid memory leaks. If you want to run in a single thread, e.g. to be able to read the debug output during generation, you should set the option to zero.
If running multi-threaded on the cern batch cluster, you should declare it to the `bsub` command when submitting the jobs: e.g. for a job that uses four cores you should use **`bsub -n 4 -R "'span[hosts=1]'" ...`**

\#HybridNewGrid

##### HybridNew algorithm usage for complex models or expected limits: grids

If your model is complex, or you need to know the limit accurately, or you want expected limits, then running the computation in a single job might not be feasible.

The alternative approach is to compute a grid of distributions of the test statistics for various values of the signal strength, a task that is easy to parallelize, and then use the that grid to compute the observed limit (and also the expected ones). This requires you to have some knowledge of where the limit should be, which you can gain e.g. from the ProfileLikelihood method

**Creating the grid: manual way**
The procedure to do this manually would be like the procedure for significances or simple hypothesis testing described previously: for each value **`r_i`** of the cross section, you write out one file with the distribution of the test statistics using

    combine card.txt -M HybridNew [--freq] [other options] -s seed_i --singlePoint r_i --saveToys --saveHybridResult

and then you can merge all the output files for the different **`r_i`** with **`hadd`**. The `[other options]` should include **`--clsAcc 0`** to switch off adaptive sampling, and you can tune the CPU time by working on the parameters **`-T`** and **`--iterations`**.
It is important that you use different **`seed_i`** values for each point; if you don't care about exact reproducibility, you can just use **`--seed -1`** and the code will take care of randomizing itself properly.

**Creating the grid: automated way, using CRAB**

Please note that the following is intended for use with crab2. For producing the grid with crab3, please see the instructions [here](https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/SWGuideNonStandardCombineUses#Running_combine_on_the_grid)

Once you have a sense of the time needed for each toy, and of the range to consider, you can use the script [makeGridUsingCrab.py](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/test/makeGridUsingCrab.py) to run the toys to create the grid in parallel either on LXBATCH or on regular T2s (or anything else that CRAB can digest). The procedure is very simple:

    makeGridUsingCrab.py card.txt minimum maximum -n points [other options] -o name

This will create a crab cfg file `name.cfg` and a script `name.sh`, and possibly a binary workspace `name.workspace.root`. You can then just create and submit the jobs from that cfg file, and merge the output rootfiles with `hadd=(note:  =hadd` will complain with messages like `Cannot merge object type, name: HybridCalculator _result` which can be safely ignored).
The other options, that you can get executing the program with **`--help`** are:

-   **`-T`**: same as in combine
-   **`-r`**: use a random seed in each job (suggested)
-   **`-I n`**: run only on 1/n of the points in each job (suggested if you want to have many points)
-   **`-t`**, **`-j`**: choose the total number of toys and of jobs (can change later from the crab cfg file)
-   **`--lsf`**, **`--queue ...`**: use lxbatch with the specific queue (can change later from the crab cfg file)

Note that you can merge also the output of multiple crab submissions, if you have used random seeds.

**Using the grid for observed limits**

    combine mydatcard.txt -M HybridNew [--freq] --grid=mygrid.root

All other parameters controlling toys, accuracy and so on are meaningless in this context. Note that it might still take a while if you have many points and the test statistics is slow to evaluate. Add the option **`--saveGrid`** to save the value of the observed CLs at each grid point in the output tree.

**Using the grid for expected limits**

    combine mydatcard.txt -M HybridNew [--freq] --grid=mygrid.root --expectedFromGrid 0.5

0.5 gives you the median. use 0.16/0.84 to get the endpoints of 68% interval, 0.025/0.975 to get the 95% one). Add the option **`--saveGrid`** to save the value of the expected quantile CLs at each grid point in the output tree.

**Plotting the test-statistics distributions**

The distribution of the test-statistic under the signal plus background and background only hypotheses can be plotted at each value of the grid using the following;

    python test/plotTestStatCLs.py --input mygrid.root --poi r --val all --mass MASS

The output root file will contain a plot for each point found in the grid.

#### FeldmanCousins

The F-C method is used to compute an interval with the specified confidence level.

If you run the model without special options, it will report the upper limit to the signal strength. If you want instead the lower end of the interval, just run it with option **`lowerLimit`**.

The algorithm will search for a limit with an iterative procedure until the specified absolute or relative accuracy is met, as controlled by the parameters **`rAbsAcc`**, `=rRelAcc`. The default values are 0.1 and 0.05 respectively, meaning that the search is stopped when Δr &lt; 0.1 or Δr/r &lt; 0.05.

The number of toys generated is adaptive. You can increase it by a factor using option **`toysFactor`** a value of 2 or 4 is suggested if you want to compute the limit with high accuracy.

Running under CRAB
------------------

<span class="twiki-macro RED"></span> The instructions below are for use with **`crab2`**. For instructions on how to use the grid for toy based studies or complicated model scans under **`crab3`**, follow the instructions given [here](https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/SWGuideNonStandardCombineUses#Running_combine_on_the_grid).

Running many toy MC for the limit calculation may be conveniently split among the different available GRID resources using CRAB. Examples of how to run on the GRID via CRAB are provided in the files:

-   **`[[https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/test/combine_crab.sh][combine_crab.sh]]`**
-   **`[[https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/test/combine_crab.cfg][combine_crab.cfg]]`**

### Preparing the ROOT workspace

The first thing to do is to convert the datacards and possibly the shape model into a ROOT workspace. This model will be shipped to Worker Nodes for GRID processing. This is done via the utility script **`text2workspace.py`**. For instance:

    ../scripts/text2workspace.py ../data/benchmarks/simple-counting/counting-B0-Obs0-StatOnly.txt -b -o model.root

### Shell script for GRID Worker Nodes

CRAB is designed mainly to provide automatic **`cmsRun`** job splitting providing the number of jobs one wants to run, and the number of 'events' in total one wants to process.The total number of toy MC we want to run. The maximum number of events is passed to the application to be executed via the variable **`$MaxEvents`**. In our case, we will use for convenience **`$MaxEvents`** as the number of toy MC to run per job.

The script **`[[https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/test/combine_crab.sh][combine_crab.sh]]`** runs the combiner code with the proper options, and prepares the output to be retrieved after the run completion on the Worker Nodes. It takes as argument the job indes (**`$1`**), which we use as random seed. The main elements there are running the combiner and packing the output for final retrieval: <span class="twiki-macro SYNTAX" syntax="sh"></span> echo "job number: seed \#$i with $n toys" ./combine model.root -t$n -s$i &gt;& log.txt mv \*.root outputToy/ mv log.txt outputToy/ echo "pack the results" tar cvfz outputToy.tgz outputToy/ <span class="twiki-macro ENDSYNTAX"></span>

### CRAB configuration script

Then, the script `= [[https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/master/test/combine_crab.cfg][combine_crab.cfg]]=` sets how many jobs to run and how many toys per job. And finally, which files to ship to the Worker Nodes (the executable itself and the ROOT workspace), and which file to retrieve (the packed output): <span class="twiki-macro SYNTAX" syntax="sh"></span> \[CRAB\] jobtype = cmssw scheduler = glite

\[CMSSW\] output\_file = outputToy.tgz datasetpath=None pset=None total\_number\_of\_events=100 number\_of\_jobs=10

\[USER\] script\_exe = combine\_crab.sh additional\_input\_files = combine, model.root return\_data = 1 <span class="twiki-macro ENDSYNTAX"></span>

Please, notice that you have to ship to Worker Nodes the executable itself, so, before running you have to do:

    ln -s ../../../../bin/slc5_amd64_gcc434/combine combine

### GRID submission via CRAB

CRAB submission can be now performed as usually done for other CMSSW analysis application:

-   create your jobs:

<!-- -->

    crab -create -cfg combine_crab.cfg

-   submit your jobs:

<!-- -->

    crab -submit

-   monitor your jobs' status:

<!-- -->

    crab -status

-   after jobs are finished, retrieve the output:

<!-- -->

    crab -getoutput

The output consists in several files `outputToy_n_m_xyz.tgz` that contain the output of each job, including the log files, that can be combined and analyzed to obtain the final result.

Useful Links and tutorials
--------------------------

-   Tutorial Sessions
    -   [1st tutorial 17th Nov 2015](https://indico.cern.ch/event/456547/).
    -   [2nd tutorial 30th Nov 2016](https://indico.cern.ch/event/577649/#b-229590-higgs-combine-tool-mi).

<!-- -->

-   [Worked examples from Higgs analyses](https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideCMSDataAnalysisSchool2014HiggsCombPropertiesExercise)

<!-- -->

-   [Further Reading and Advanced usage tutorials](CMS/HiggsWG/SWGuideNonStandardCombineUses)

<!-- -->

-   [Conventions to be used when preparing inputs for Higgs combinations](CMS/HiggsWG/HiggsCombinationConventions)

<!-- -->

-   [CMS AN-2011/298: Procedure for the LHC Higgs boson search combination in summer 2011](http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2011/298) (describing some of the methods used in `combine`)

Combine-based packages
----------------------

-   SWGuideHiggs2TauLimits

<!-- -->

-   ATGCRooStats

<!-- -->

-   [CombineHarvester](http://cms-analysis.github.io/CombineHarvester/)

Contacts
--------

-   **Hypernews forum**: hn-cms-higgs-combination <https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination.html>

\#ReviewStatus

Review status
-------------

| Reviewer/Editor and Date (copy from screen) | Comments                              |
|:--------------------------------------------|:--------------------------------------|
| Main.GiovanniPetrucciani - 14-Dec-2010      | created template page                 |
| Main.GiovanniPetrucciani - 21-Mar-2011      | new tag                               |
| Main.LucaLista - 04-May-2011                | added CRAB submission                 |
| Main.GiovanniPetrucciani - 10-Jun-2011      | Frequentist limits                    |
| Main.GiovanniPetrucciani - 10-Jun-2011      | Frequentist grids and expected limits |

<span class="twiki-macro RESPONSIBLE"></span> Main.GiovanniPetrucciani
<span class="twiki-macro REVIEW"></span> Main.GiovanniPetrucciani
