# Groups of nuisances, snapshots, and uncertainty breakdown by freezing

The files in this directory are a tutorial on:
- defining groups of nuisances,
- creating snapshots from fits, and
- breaking down uncertainties by freezing groups of nuisances while starting from the snaphot.

The files in this directory are:
- [`myanalysis.dc.txt`](./myanalysis.dc.txt): a two-bin counting experiment datacard with some tension built into the observation and multiple nuisance groups defined.
- [`Makefile`](./Makefile): an annotated list of all the commands used, as well as the walkthrough of the results.
- [`MLFUncDiff.py`](./MLFUncDiff.py): a (rather dump) python script that takes two output files from `MaxLikelihoodFit` and subtracts uncertainties in quadrature.

For simplicity, this tutorial does not discuss scans using `-M MultiDimFit --algo=grid`,
but rather uses `-M MaxLikelihoodFit --robustFit=1`
to determine the crossings of the likelihood with the 1 sigma levels.

## Groups of nuisances

Groups of nuisances are defined in the datacard, e.g. [here](./myanalysis.dc.txt#L30), [here](./myanalysis.dc.txt#L38-L40), [here](./myanalysis.dc.txt#L46-L47), and [here](./myanalysis.dc.txt#L56-L58).
They are parsed when building the binary datacard and result in two actions:
- creation in the workspace of a `RooArgSet` for each group, containing the concerned nuisance parameters.
- setting of a `group_%s` attribute in each nuisance parameter.

A nuisance parameter can belong to multiple groups of nuisances.
It is up to you to list or not list nuisances as part of groups in an exclusive, or overlapping way.
In particular, there is no check that every nuisance is in at least one group.

At present, the syntax does not allow to declare groups as members of other groups;
only individual nuisances can be declared as part of a group.

## Snapshots

Snapshots are a functionality of `MultiDimFit` that allows to save the result of the fit, see [here](./Makefile#L115-L118).
They are essential to extract post-fit uncertainties,
since the nuisances have to be set to their best-fit values.

For instance, when determining the statistical component of the uncertainty,
all nuisance parameters have to be frozen.
If nuisance parameters were frozen to their default values,
the knowledge gained on the nuisance parameters
by fitting the data would not be taken into account and even the central value would be different.

## Freezing nuisances and breaking down the uncertainty

One way to deconstruct the total uncertainty into different components is by freezing nuisances.
For instance, freezing all nuisances directly gives the statistical component as mentioned above.
By subtracting in quadrature the statistical uncertainty from the total uncertainty, one can derive the
systematic uncertainty.

When there are multiple components being disentangled,
the decomposition can proceed in many ways,
and an example of extracting the exp. syst. and theo. syst. is given in this tutorial, [here](./Makefile#L42-L70)

## Show me the money

To run the tutorial, you need a properly set up `combine` area and run `make` in this directory.
The output walks you through the whole procedure and should look something like this:
```ShellSession
$ make
Creating binary datacard:
text2workspace.py -o myanalysis.ws.root -v9 myanalysis.dc.txt &> logs/myanalysis.T2W.log
Running and saving the best-fit snapshot:
combine -M MultiDimFit --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 --saveWorkspace -n _myanalysis_bestfit -v9 myanalysis.ws.root &> logs/myanalysis_bestfit.log
Profiling with [bgexp] group of nuisances frozen to best-fit values:
combine -M MaxLikelihoodFit --freezeNuisanceGroups bgexp        --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 -d higgsCombine_myanalysis_bestfit.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _frozen_bgexp -v9 &> logs/frozen_bgexp.log
Profiling with [exp] group of nuisances frozen to best-fit values:
combine -M MaxLikelihoodFit --freezeNuisanceGroups exp        --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 -d higgsCombine_myanalysis_bestfit.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _frozen_exp -v9 &> logs/frozen_exp.log
Profiling with [sigtheo] group of nuisances frozen to best-fit values:
combine -M MaxLikelihoodFit --freezeNuisanceGroups sigtheo        --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 -d higgsCombine_myanalysis_bestfit.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _frozen_sigtheo -v9 &> logs/frozen_sigtheo.log
Profiling with [bgtheo] group of nuisances frozen to best-fit values:
combine -M MaxLikelihoodFit --freezeNuisanceGroups bgtheo        --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 -d higgsCombine_myanalysis_bestfit.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _frozen_bgtheo -v9 &> logs/frozen_bgtheo.log
Profiling with [theo] group of nuisances frozen to best-fit values:
combine -M MaxLikelihoodFit --freezeNuisanceGroups theo        --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 -d higgsCombine_myanalysis_bestfit.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _frozen_theo -v9 &> logs/frozen_theo.log
Profiling with [bg] group of nuisances frozen to best-fit values:
combine -M MaxLikelihoodFit --freezeNuisanceGroups bg        --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 -d higgsCombine_myanalysis_bestfit.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _frozen_bg -v9 &> logs/frozen_bg.log
Profiling with all nuisances frozen to best-fit values:
#TODO --profilingMode none is not working?
combine -M MaxLikelihoodFit --freezeNuisanceGroups exp,theo  --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 -d higgsCombine_myanalysis_bestfit.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _frozen_all -v9 &> logs/frozen_all.log
Profiling with no nuisances frozen:
combine -M MaxLikelihoodFit                                  --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 -d higgsCombine_myanalysis_bestfit.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _frozen_none -v9 &> logs/frozen_none.log
Profiling with all nuisances frozen (and not using the snapshot, so frozen to default values):
#TODO --profilingMode none is not working?
combine -M MaxLikelihoodFit --freezeNuisanceGroups exp,theo  --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 myanalysis.ws.root                                      -n _frozen_all_mlf -v9 &> logs/frozen_all_mlf.log
Profiling with no nuisances frozen (and not using the snapshot):
combine -M MaxLikelihoodFit                                  --robustFit 1 --minimizerTolerance 0.00000001 --minimizerStrategy 2 myanalysis.ws.root                                      -n _frozen_none_mlf -v9 &> logs/frozen_none_mlf.log

Welcome to a weird way of presenting results that is halfway between slides and jupyter.

0) Note the text2workspace printout regarding groups, that helps figuring out what is going on:
group_bg:(ZZnorm,ZHnorm,DY_NLO,tWnorm,bkgNorm2,bkgNorm3,bkgNorm1,WZnorm,WWnorm)
group_bgexp:(bkgNorm2,bkgNorm3,bkgNorm1)
group_bgtheo:(ZZnorm,ZHnorm,WZnorm,tWnorm,DY_NLO,WWnorm)
group_exp:(CMS_scale_j,CMS_res_j,El_SFs,lumi_8TeV,Mu_SFs,btagSFs_bcjets,bkgNorm3,bkgNorm1,bkgNorm2,btagSFs_lightjets)
group_sigtheo:(SIGpfd,SIGeff,SIGscale)
group_theo:(ZZnorm,ZHnorm,DY_NLO,SIGpfd,SIGeff,SIGscale,tWnorm,WZnorm,WWnorm)

1) The best-fit snaphot has the following best-fit values:
CMS_res_j	  = -0.132802	 +/-  0.989801
CMS_scale_j	  = -0.173259	 +/-  0.976373
DY_NLO	  = 0.0191775	 +/-  0.994808
El_SFs	  = -0.392716	 +/-  0.921826
Mu_SFs	  = 0.412567	 +/-  0.913366
SIGeff	  = 1.34081e-08	 +/-  0.994458
SIGpfd	  = 3.95133e-08	 +/-  0.994465
SIGscale	  = 3.39636e-08	 +/-  0.994463
WWnorm	  = 0.000241014	 +/-  0.994468
WZnorm	  = 0.000685482	 +/-  0.994476
ZHnorm	  = 0.000440988	 +/-  0.994471
ZZnorm	  = 0.00281734	 +/-  0.994601
bkgNorm1	  = -0.0192619	 +/-  0.993785
bkgNorm2	  = -0.00167932	 +/-  0.994539
bkgNorm3	  = 5.724e-06	 +/-  0.994672
bkgNorm4	  = -0.0111463	 +/-  0.994518
btagSFs_bcjets	  = -0.00518017	 +/-  0.994563
btagSFs_lightjets	  = 0.0343742	 +/-  0.995218
lumi_8TeV	  = 0.00274212	 +/-  0.994492
r	  = 0.935189	 +/-  0.631912	(limited)
tWnorm	  = 0.0121632	 +/-  0.995646
Minimized in : Real time 0:00:02, CP time 0.550
RooRealVar::r = 0.935189 +/- 0.631912  L(0 - 20)

 --- MultiDimFit ---
best fit parameter values:
   r :    +0.935
Done in 0.01 min (cpu), 0.04 min (real)

=> Note how nuisances are pulled from zero.

2) Now, compare what happens if you freeze nuisances to the default values directly in MLF (all_mlf) or to the best-fit values using the snapshot (all):
==> logs/frozen_all_mlf.log <==
Best fit r: 0.956118  -0.40413/+0.409622  (68% CL)
--
==> logs/frozen_all.log <==
Best fit r: 0.935188  -0.406103/+0.412308  (68% CL)

=> It makes a big difference!

3) This difference does not exist for the unconstrained fit in MLF (none_mlf) or using the snapshot (none):
==> logs/frozen_none_mlf.log <==
Best fit r: 0.935189  -0.628684/+0.646119  (68% CL)
--
==> logs/frozen_none.log <==
Best fit r: 0.935186  -0.628684/+0.646123  (68% CL)

=> This is just a (successful) closure test.

4) What is the systematic uncertainty?
   Assuming (in quadrature notation) that tot = stat + syst, we can estimate syst = tot - stat:

4.1) The (post-fit) full uncertainty is given by (frozen_none):
Best fit r: 0.935186  -0.628684/+0.646123  (68% CL)

4.2) The (post-fit) stat-only uncertainty is given by (frozen_all):
Best fit r: 0.935188  -0.406103/+0.412308  (68% CL)

4.3) The difference (the syst) is:
logs/frozen_none.log:Best fit r: 0.935186  -0.628684/+0.646123  (68% CL)
logs/frozen_all.log:Best fit r: 0.935188  -0.406103/+0.412308  (68% CL)
(-0.4799202495512123, 0.49747090428258917)

5) Let's get the exp. syst. uncertainty.
   Assuming (in quadrature notation) that tot = stat + exp + theo, and given the way we defined the groups of nuisances, there are two ways:
   - (exp) = (tot) - (stat+theo), or
   - (exp) = (tot-theo) - (stat).
  From the fits, we have all the terms above in brackets:

5.1) (exp) = (tot) - (stat+theo), which translates to (exp) = [frozen_none] - [frozen_exp]
logs/frozen_none.log:Best fit r: 0.935186  -0.628684/+0.646123  (68% CL)
logs/frozen_exp.log:Best fit r: 0.935187  -0.555647/+0.560975  (68% CL)
(-0.2941083488468959, 0.32059721177096867)

5.2) (exp) = (tot-theo) - (stat), which translates to (exp) = [frozen_theo] - [frozen_all]
logs/frozen_theo.log:Best fit r: 0.935186  -0.50146/+0.520246  (68% CL)
logs/frozen_all.log:Best fit r: 0.935188  -0.406103/+0.412308  (68% CL)
(-0.29418057987173946, 0.3172674850846754)

=> These are very similar as expected.

6) Let's get the theo. syst. uncertainty.
   Assuming (in quadrature notation) that tot = stat + exp + theo, and given the way we defined the groups of nuisances, there are two ways:
   - (theo) = (tot) - (stat+exp), or
   - (theo) = (tot-exp) - (stat).
  From the fits, we have all the terms above in brackets:

6.1) (theo) = (tot) - (stat+exp), which translates to (theo) = [frozen_none] - [frozen_theo]
logs/frozen_none.log:Best fit r: 0.935186  -0.628684/+0.646123  (68% CL)
logs/frozen_theo.log:Best fit r: 0.935186  -0.50146/+0.520246  (68% CL)
(-0.37918495797648016, 0.3831692100310022)

6.2) (theo) = (tot-exp) - (stat), which translates to (theo) = [frozen_exp] - [frozen_all]
logs/frozen_exp.log:Best fit r: 0.935187  -0.555647/+0.560975  (68% CL)
logs/frozen_all.log:Best fit r: 0.935188  -0.406103/+0.412308  (68% CL)
(-0.3792409854800118, 0.3803876028637338)

=> These are also very similar as expected.

7) Let's look at partial effects of groups of nuisances.
   For this case we will proceed with (group) = (tot) - (stat+group)

7.1)
group_bg:(ZZnorm,ZHnorm,DY_NLO,tWnorm,bkgNorm2,bkgNorm3,bkgNorm1,WZnorm,WWnorm)
logs/frozen_none.log:Best fit r: 0.935186  -0.628684/+0.646123  (68% CL)
logs/frozen_bg.log:Best fit r: 0.935187  -0.487043/+0.517858  (68% CL)
(-0.3975327938103605, 0.38639159547611013)

7.2)
group_bgtheo:(ZZnorm,ZHnorm,WZnorm,tWnorm,DY_NLO,WWnorm)
logs/frozen_none.log:Best fit r: 0.935186  -0.628684/+0.646123  (68% CL)
logs/frozen_bgtheo.log:Best fit r: 0.935186  -0.502894/+0.532303  (68% CL)
(-0.37728070424293597, 0.3662354114117034)

7.2)
group_bgexp:(bkgNorm2,bkgNorm3,bkgNorm1)
logs/frozen_none.log:Best fit r: 0.935186  -0.628684/+0.646123  (68% CL)
logs/frozen_bgexp.log:Best fit r: 0.935186  -0.616016/+0.634923  (68% CL)
(-0.12557002567380163, 0.11978438405015161)

7.4)
group_sigtheo:(SIGpfd,SIGeff,SIGscale)
logs/frozen_none.log:Best fit r: 0.935186  -0.628684/+0.646123  (68% CL)
logs/frozen_sigtheo.log:Best fit r: 0.935186  -0.628185/+0.634439  (68% CL)
(-0.02503213081439212, 0.12231861977354595)

8) What can we learn?
   With these, we see that this analysis would be syst-limited (4.3 vs 4.2),
   with the most important systematics being from theory (6 vs 5),
   namely the theory part on the bg (7.2 vs 5 and 7.3).
   To learn more you need to look at the individual impacts <insert here next lecture>.
```
