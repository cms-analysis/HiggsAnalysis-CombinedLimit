# Answers to tasks and questions in long exercise

## Part 1: A one-bin counting experiment

### A: Computing limits using the asymptotic approximation

**Tasks and questions:**

  -   There are some important uncertainties missing from the datacard above. Add the uncertainty on the luminosity (name: `lumi_13TeV`) which has a 2.5% effect on all processes (except the `jetFakes`, which are taken from data), and uncertainties on the inclusive cross sections of the `Ztautau` and `ttbar` processes (with names `xsec_Ztautau` and `xsec_diboson`) which are 4% and 6% respectively.
  -   Try changing the values of some uncertainties (up or down, or removing them altogether) - how do the expected and observed limits change?
<details>
<summary><b>Show answer</b></summary>

<b> Larger uncertainties make the limits worse (ie, higher values of the limit); smaller uncertainties improve the limit (lower values of the limit). </b>

</details>

  -   Now try changing the number of observed events. The observed limit will naturally change, but the expected does too - why might this be?
<details>
<summary><b>Show answer</b></summary>

<b> This is because the expected limit relies on a background-only Asimov dataset that is created <em> after </em> a background-only fit to the data. By changing the observed the pulls on the NPs in this fit also change, and therefore so does the expected sensitivity.</b>

</details>

### Advanced section: B: Computing limits with toys

**Tasks and questions:**

  - In contrast to `AsymptoticLimits`, with `HybridNew` each limit comes with an uncertainty. What is the origin of this uncertainty?
<details>
<summary><b>Show answer</b></summary>

<b> The uncertainty is caused by the limited number of toys: the values of Pmu and Pb come from counting the number of toys in the tails of the test statistic distributions. The number of toys used can be adjusted with the option <code> --toysH </code> </b>

</details>
  - How good is the agreement between the asymptotic and toy-based methods?
<details>
<summary><b>Show answer</b></summary>

<b> The agreement should be pretty good in this example, but will generally break down once we get to the level of 0-5 events. </b>

</details>
  - Why does it take longer to calculate the lower expected quantiles (e.g. 0.025, 0.16)? Think about how the statistical uncertainty on the CLs value depends on Pmu and Pb.
<details>
<summary><b>Show answer</b></summary>

<b> For this we need the definition of CLs = Pmu / (1-Pb). The 0.025 expected quantile is by definition where Pb = 0.025, so for a 95% CL limit we have CLs = 0.05, implying we are looking for the value of r where Pmu = 0.00125. With 1000 s+b toys we would then only expect `1000 * 0.00125 = 1.25` toys in the tail region we have to integrate over. Contrast this to the median limit where 25 toys would be in this region. This means we have to generate a much larger numbers of toys to get the same statistical power.</b>

</details>

### Advanced section: B: Asymptotic approximation limitations

**Tasks and questions:**

  - Is the asymptotic limit still a good approximation?
<details>
<summary><b>Show answer</b></summary>

  <b> A "good" approximation is not well defined, but the difference is clearly larger here. </b>

</details>
  - You might notice that the test statistic distributions are not smooth but rather have several "bump" structures? Where might this come from? Try reducing the size of the systematic uncertainties to make them more pronounced.
<details>
<summary><b>Show answer</b></summary>

<b> This bump structure comes from the discrete-ness of the Poisson sampling of the toy datasets. Systematic uncertainties then smear these bumps out, but without systematics we would see delta functions corresponding to the possible integer number of events that could be observed. Once we go to more typical multi-bin analyses with more events and systematic uncertainties these discrete-ness washes out very quickly. </b>

</details>


## Part 2: A shape-based analysis

### A: Setting up the datacard

Only tasks, no questions in this section

### B: Running combine for a blind analysis

**Tasks and questions:**

  - Compare the expected limits calculated with --run expected and --run blind. Why are they different?
<details>
<summary><b>Show answer</b></summary>

<b> When using --run blind combine will create a background-only Asimov dataset without performing a fit to data first. With --run expected, the observed limit isn't shown, but the background-only Asimov dataset used for the limit calculation is still created after a background-only fit to the data. </b>

</details>
  - Calculate a blind limit by generating a background-only Asimov with the -t option instead of using the AsymptoticLimits specific options. You should find the observed limit is the same as the expected. Then see what happens if you inject a signal into the Asimov dataset using the --expectSignal [X] option.
<details>
<summary><b>Show answer</b></summary>

<b> You should see that with a signal injected the observed limit is worse (has a higher value) than the expected limit: for the expected limit the b-only Asimov dataset is still used, but the observed limit is now calculated on the signal + background Asimov dataset, with a signal at the specified cross section [X]. </b>

</details>

### C: Using FitDiagnostics

**Tasks and questions:**

  - Which parameter has the largest shift from the nominal value? Which has the tightest constraint?
<details>
<summary><b>Show answer</b></summary>
  
<b>  <code> CMS_eff_t_highpt </code> should have the largest shift from the nominal value (around 0.47), <code> norm_jetFakes </code> has the tightest constraint (to 25% of the input uncertainty). </b>

</details>
  - Should we be concerned when a parameter is more strongly constrained than the input uncertainty (i.e. $\frac{\sigma}{\sigma_I}<1.0$)?
<details>
<summary><b>Show answer</b></summary>

<b> This is still a hot topic in CMS analyses today, and there isn't a right or wrong answer. Essentially we have to judge if our analysis should really be able to provide more information about this parameter than the external measurement that gave us the input uncertainty. So we would not expect to be able to constrain the luminosity uncertainty for example, but uncertainties specific to the analysis might legitimately be constrained. </b>

</details>


### D: MC statistical uncertainties
**Tasks and questions:**

  - Check how much the cross section measurement and uncertainties change using `FitDiagnostics`.
<details>
<summary><b>Show answer</b></summary>

  <b> Without autoMCStats we find: <code> Best fit r: -2.73273  -2.13428/+3.38185</code>, with autoMCStats: <code> Best fit r: -3.07825  -3.17742/+3.7087 </code> </b>

</details>
  - It is also useful to check how the expected uncertainty changes using an Asimov dataset, say with `r=10` injected.
<details>
<summary><b>Show answer</b></summary>

  <b> Without autoMCStats we find: <code> Best fit r: 9.99978  -4.85341/+6.56233 </code>, with autoMCStats: <code> Best fit r: 9.99985  -5.24634/+6.98266 </code> </b>

</details>
  - **Advanced task:** See what happens if the Poisson threshold is increased. Based on your results, what threshold would you recommend for this analysis?
<details>
<summary><b>Show answer</b></summary>

<b> At first the uncertainties increase, as the threshold increases, and at some point they stabilise. A Poisson threshold at 10 is probably reasonable for this analysis. </b>

</details>

## Part 3: Adding control regions

### A: Use of rateParams

**Tasks and questions:**

  - Run `text2workspace.py` on this combined card and then use `FitDiagnostics` on an Asimov dataset with `r=1` to get the expected uncertainty. Suggested command line options: `--rMin 0 --rMax 2`
<details>
<summary><b>Show answer</b></summary>

  <b> As expected uncertainty you should get <code> -0.417238/+0.450593 </code> </b>

</details>
  - Using the RooFitResult in the `fitDiagnosticsTest.root` file, check the post-fit value of the rateParams. To what level are the normalisations of the DY and ttbar processes constrained?
<details>
<summary><b>Show answer</b></summary>

  <b> They are constrained to around 1-2% </b>

</details>
  - To compare to the previous approach of fitting the SR only, with cross section and acceptance uncertainties restored, an additional card is provided: `datacard_part3_nocrs.txt`. Run the same fit on this card to verify the improvement of the SR+CR approach
<details>
<summary><b>Show answer</b></summary>

  <b> The expected uncertainty is larger with only the SR: <code> -0.465799/+0.502088 </code> compared with <code> -0.417238/+0.450593 </code> in the SR+CR approach. </b>

</details>


### B: Nuisance parameter impacts
**Tasks and questions:**

  - Identify the most important uncertainties using the impacts tool.
<details>
<summary><b>Show answer</b></summary>
  
 <b> The most important uncertainty is <code>norm_jetFakes</code>, followed by two MC statistical uncerainties (<code>prop_binsignal_region_bin8</code> and <code>prop_binsignal_region_bin9</code>).</b>

</details>
  - In the plot, some parameters do not show a plotted point for the fitted value, but rather just a numerical value - why?
<details>
<summary><b>Show answer</b></summary>

  <b> These are freely floating parameters (<code> rate_ttbar </code> and <code> rate_Zll </code>). They have no prior constraint (and so no shift from the nominal value relative to the input uncertainty) - we show the best-fit value + uncertainty directly. </b>

</details>


### C: Post-fit distributions
**Tasks and questions:**

 The bin errors on the TH1s in the fitdiagnostics file are determined from the systematic uncertainties. In the post-fit these take into account the additional constraints on the nuisance parameters as well as any correlations.

  - Why is the uncertainty on the post-fit so much smaller than on the pre-fit?
<details>
<summary><b>Show answer</b></summary>

<b> There are two effects at play here: the nuisance parameters get constrained, and there are anti-correlations between the parameters which also have the effect of reducing the total uncertainty. Note: the post-fit uncertainty could become larger when rateParams are present as they are not taken into account in the pre-fit uncertainty but do enter in the post-fit uncertainty. </b>

</details>


### D: Calculating the significance
**Tasks and questions:**

  - **Advanced task** It is also possible to calculate the significance using toys with `HybridNew` (details [here](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/#computing-significances-with-toys)) if we are in a situation where the asymptotic approximation is not reliable or if we just want to verify the result. Why might this be challenging for a high significance, say larger than $5\sigma$?
<details>
<summary><b>Show answer</b></summary>

<b> A significance of $5\sigma$ corresponds to a p-value of around $3\cdot 10^{-7}$ - so we need to populate the very tail of the test statistic distribution and this requires generating a large number of toys. </b>

</details>


### E: Signal strength measurement and uncertainty breakdown

** Tasks and questions: **

  - Take our stat+syst split one step further and separate the systematic part into two: one part for hadronic tau uncertainties and one for all others. Do this by defining a `tauID` group in the datacard including the following parameters: `CMS_eff_t`, `CMS_eff_t_highpt`, and the three `CMS_scale_t_X` uncertainties.

<details>
<summary><b>Show datacard line</b></summary>

<b> You should add this line to the end of the datacard:
```shell
tauID group = CMS_eff_t CMS_eff_t_highpt CMS_scale_t_1prong0pi0_13TeV CMS_scale_t_1prong1pi0_13TeV CMS_scale_t_3prong0pi0_13TeV
```
  </b>
</details>
  - To plot this and calculate the split via the relations above you can just add further arguments to the `--others` option in the `plot1DScan.py` script. Each is of the form: `'[file]:[label]:[color]'`. The `--breakdown` argument should also be extended to three terms.

<details>
<summary><b>Show code</b></summary>

<b> This can be done as:
```shell
python plot1DScan.py higgsCombine.part3E.MultiDimFit.mH200.root --others 'higgsCombine.part3E.freezeTauID.MultiDimFit.mH200.root:FreezeTauID:4' 'higgsCombine.part3E.freezeAll.MultiDimFit.mH200.root:FreezeAll:2' -o freeze_third_attempt --breakdown TauID,OtherSyst,Stat

```
  </b>
</details>
  - How important are these tau-related uncertainties compared to the others?
<details>
<summary><b>Show answer</b></summary>

<b> They are smaller than both the statistical uncertainty and the remaining systematic uncertainties </b>

</details>


### F: Use of channel masking

No specific questions, just tasks
