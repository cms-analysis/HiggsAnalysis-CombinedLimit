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

## Part 1: Setting up the datacard and the workspace

Topics covered in this section:

  - A: Setting up the datacard and the workspace
  - B: MC statistical uncertainties
  - C: Running <span style="font-variant:small-caps;">Combine</span> for a blind analysis -> TO ADD: run with/without shapes and compare sensitivity
  - D: Using FitDiagnostics to validate your setup
  - Extra: CAT gitLab tools for validation

### A: Setting up the datacard
In a typical analysis we will produce some distribution of our observables, so that we can separate signal and background processes and compare them with data. Combine can receive as input TH1, RooHists, and RooPDFs of any dimensions. In this example we will use TH1 histograms: one for the data and one for each signal and background processes. An example datacard should look like:

<details>
<summary><b>Show datacard</b></summary>
```shell
imax 1
jmax 1
kmax *
---------------
shapes * * simple-shapes-TH1_input.root $PROCESS $PROCESS_$SYSTEMATIC
shapes signal * simple-shapes-TH1_input.root $PROCESS$MASS $PROCESS$MASS_$SYSTEMATIC
---------------
bin bin1
observation 85
------------------------------
bin             bin1       bin1
process         signal     background
process         0          1
rate            10         100
--------------------------------
lumi     lnN    1.10       1.0
bgnorm   lnN    1.00       1.3
alpha  shape    -          1
```
</details>

The first block tells Combine (and readers) the number of bins/observables (imax), the number of background processes (jmax) and the number of nuisance parameters (kmax).
The second block tells Combine where it can find the input shapes, according to the pattern `shapes [process] [channel] [file] [histogram] [histogram_with_systematics]`. It is possible to use the `*` wildcard to map multiple processes and/or channels with one line. The histogram entries can contain the `$PROCESS`, `$CHANNEL` and `$MASS` place-holders which will be substituted when searching for a given (process, channel) combination. The value of `$MASS` is specified by the `-m` argument when combine. The final argument of the "shapes" line above should contain the `$SYSTEMATIC` place-holder which will be substituted by the systematic name given in the datacard. By default the observed data process name will be `data_obs`.
The third block lists for each bin the processes contributiong to it and their rates (yields)
Finally, the last block lists the nuisance parameters. a lnN uncertainty means that the nuisance only affects the overall normalization, while a `shape` uncertainty affects the distributions of the events.

Note that the total nominal rate of a given process, specified in the `rate` line of the datacar, must agree with the value returned by `TH1::Integral`. However, we can also put a value of `-1` and the Integral value will be substituted automatically. While it makes easier to write down the cards, this features makes debugging and reading the card harder, use it at your own risk!

Shape uncertainties can be added by supplying two additional histograms for a process, corresponding to the distribution obtained by shifting that parameter up and down by one standard deviation. These shapes will be interpolated (see the [template shape uncertainties](../part2/settinguptheanalysis.md#template-shape-uncertainties) section for details) for shifts within $\pm1\sigma$ and linearly extrapolated beyond. The normalizations are interpolated linearly in log scale just like we do for log-normal uncertainties.

![](images/shape_morphing.jpg)

In the list of uncertainties the interpretation of the values for `shape` lines is a bit different from `lnN`. The effect can be "-" or 0 for no effect, 1 for normal effect, and possibly something different from 1 to test larger or smaller effects (in that case, the unit Gaussian is scaled by that factor before using it as parameter for the interpolation).


CAREFUL FROM HERE IS EXERCISE! I NEED TO VERIFY THE FLOW. Maybe we can even skip this part and provide directly the card?

**Tasks and questions:**

Have a look at `datacard_part2.txt`: this is currently set up as a one-bin counting experiment, this means only yields are provided.

The **first task** is to convert this to a shape analysis: the file `datacard_part2.shapes.root` contains all the necessary histograms, including those for the relevant shape systematic uncertainties. Add the relevant `shapes` lines to the top of the datacard (after the `kmax` line) to map the processes to the correct TH1s in this file. Hint: you will need a different line for the signal process, since their naming patterns are different from the background ones.

Compared to the counting experiment we must also consider the effect of uncertainties that change the shape of the distribution. Some, like `CMS_eff_t_highpt`, are already present in the datacard, as it has both a shape and normalisation effect.

Add the following shape uncertainties: `top_pt_ttbar_shape` affecting `ttbar`,the tau energy scale uncertainties `CMS_scale_t_1prong0pi0_13TeV`, `CMS_scale_t_1prong1pi0_13TeV` and `CMS_scale_t_3prong0pi0_13TeV` affecting all processes except `jetFakes`, and `CMS_eff_t_highpt` also affecting the same processes.

Once this is done you can convert the text datacard into a RooFit workspace. If we feed the datacard directly into Combine, this step will be done internally every time we run. It is a good idea to do it explicitely especially for more complex analyses, since the conversion step can take a notable amount of time. For this we use the `text2workspace.py` command:

```shell
text2workspace.py datacard_part2.txt -m 800 -o workspace_part2.root
```
And then we can verify that our setup works properly using this as input to combine: -> SHOULD WE REMOVE THIS PIECE?
```shell
combine -M AsymptoticLimits workspace_part2.root -m 800
```

  - Verify that the sensitivity of the shape analysis is indeed improved over the counting analysis in the first part. -> TO BE MOVED LATER AS OPTIONAL EXERCISE
  - *Advanced task*: You can open the workspace ROOT file interactively and print the contents: `w->Print();`. Each process is represented by a PDF object that depends on the shape morphing nuisance parameters. From the workspace, choose a process and shape uncertainty, and make a plot overlaying the nominal shape with different values of the shape morphing nuisance parameter. You can change the value of a parameter with `w->var("X")->setVal(Y)`, and access a particular pdf with `w->pdf("Z")`. PDF objects in RooFit have a [createHistogram](https://root.cern.ch/doc/master/classRooAbsReal.html#a552a08367c964e689515f2b5c92c8bbe) method that requires the name of the observable (the variable defining the x-axis) - this is called `CMS_th1x` in combine datacards. Feel free to ask for help with this!

### B: MC Statistical uncertainties


