# Building statistical models with Combine 

## Getting started
To get started, you should have a working setup of <span style="font-variant:small-caps;">Combine</span>, please follow the instructions from the [home page](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#within-cmssw-recommended-for-cms-users). Make sure to use the latest recommended release.

After setting up <span style="font-variant:small-caps;">Combine</span>, you can access the working directory for this tutorial which contains all of the inputs and scripts needed in this excercise exercise:

```shell
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/
git checkout main 
scram b -j 8
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/tutorials/model_building_2024/
```

## Exercise outline

This tutorial focuses and extends on the model building topic, and it is not going to give a full picture on the statistical methods, which are extensively covered in the [long exercise](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/part5/longexercise/) and [statistical methods exercise](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/tutorial_stat_routines/stat_routines/). 

1) Building analysis with templates 

2) Using keywords 

3) Control regions 

4) Rate parameters 

5) Extra arguments 	

6) Physics Models

## Introduction

The most general definition for the binned model likelihood can be given as 

$$ \mathcal{L} =  \mathcal{L}_\mathrm{primary} \cdot \mathcal{L}_\mathrm{auxiliary} = \prod_{c=1}^{N_c} \prod_{b=1}^{N_b^c} \mathrm{Poiss}(n_{cb}; n^\mathrm{exp}_{cb}(\vec{\mu},\vec{\nu})) \cdot \prod_{e=1}^{N_E}  p_e(y_e ; \nu_e) $$

Where $c$ indexes the channel, $b$ indexes the histogram bin, and $e$ indexes the nuisance parameter.

The generic model of the expected event count in a given bin, $n^\mathrm{exp}_{cb}$, implemented in combine for template based analyses is given by:

$$n^\mathrm{exp}_{cb} = \mathrm{max}(0, \sum_{p} M_{cp}(\vec{\mu})N_{cp}(\nu_G, \vec{\nu}_L,\vec{\nu}_S,\vec{\nu}_{\rho})\omega_{cbp}(\vec{\nu}_S) + E_{cb}(\vec{\nu}_B) ) $$

In terms of datacard structure there are several differences with respect to the counting datacard:

  - A new block of lines at the top defining how channels and processes are mapped to the histograms (more than one line can be used)
  - In the list of systematic uncertainties we now have entries marked as `shape`

The "shapes" line has to follow the following syntax: 

```shell
shapes   <process_name>   <channel_name>   <path/to/input_shape.root>   <histograms_nominal>  <histograms_with_variations>
```

To start the hands-on for this section: 
```
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/tutorials/model_building_2024/simple_shape
```

The input shapes for the first example (`datacard-2-template-analysis.txt`) are located in the `template-analysis-datacard-input.root`, it contains the observed distribution `data_obs`, the nominal histograms for each process and systematic uncertainties templates:  

```shell
root [0] 
Attaching file template-analysis-datacard-input.root as _file0...
(TFile *) 0x556fd5b0fea0
root [1] .ls
TFile**		template-analysis-datacard-input.root	
 TFile*		template-analysis-datacard-input.root	
  KEY: TH1F	background;1	Histogram of background__x
  KEY: TH1F	background_alphaUp;1	Histogram of background__x
  KEY: TH1F	background_alphaDown;1	Histogram of background__x
  KEY: TH1F	data_obs;1	Histogram of data_obs__x
  KEY: TH1F	signal;1	Histogram of signal__x
  KEY: TH1F	signal_sigmaUp;1	Histogram of signal__x
  KEY: TH1F	signal_sigmaDown;1	Histogram of signal__x
```

To define the mapping to the systematic uncertainties templates the `$SYSTEMATIC` keyword should be used, which connects the systematic uncertainties marked as `shape` type with the input shapes. 

```shell
imax 1
jmax 1
kmax 4
# ---------------
shapes	signal	ch1	template-analysis-datacard-input.root signal signal_$SYSTEMATIC
shapes	background	ch1	template-analysis-datacard-input.root background background_$SYSTEMATIC
# ---------------
bin         ch1
observation 85
# ------------------------------
bin             ch1        ch1
process         signal     background
process         0          1
rate            24         100
# --------------------------------
lumi     lnN    1.1       1.0
bgnorm   lnN    -         1.3
alpha  shape    -          1   # uncertainty in the background template.
sigma  shape    0.5        -   # uncertainty in the signal template.
```

To simplify the shape mapping line the keywords `$PROCESS`, `$CHANNEL` can be used. The `$PROCESS` keyword is associated with the processes listed in the datacard: `[signal, background]`, it is also possible to use the `*` wildcard to map multiple processes and/or channels with one line as shown below. 

```shell
imax 1
jmax 1
kmax 4
# ---------------
shapes * * template-analysis-datacard-input.root $PROCESS $PROCESS_$SYSTEMATIC
# ---------------
bin         ch1
observation 85
# ------------------------------
bin             ch1        ch1
process         signal     background
process         0          1
rate            24         100
# --------------------------------
lumi     lnN    1.1       1.0
bgnorm   lnN    -         1.3
alpha  shape    -          1   # uncertainty in the background template.
sigma  shape    0.5        -   # uncertainty in the signal template.

```

If there are more than one category it can be useful to store the input shapes corresponding to different regions in separate `TDirectory`s and use $CHANNEL keyword as shown below: 

```shell
shapes * * <input-file.root> $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC
```

## Keywords

Go to the datacards location corresponding to this section of the tutorial:
```
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/tutorials/model_building_2024/keywords
```

The datacard can also contain `$MASS` keyword, it allows to setup a single datacard for various mass points. It will be replaced with the value passed to `-m` option when running the tool. In addition, user-defined keywords can be used. Any word in the datacard `$WORD` will be replaced by `VALUE` when including the option `--keyword-value WORD=VALUE`. This option can be repeated multiple times for multiple keywords.

```
 KEY: TH1D	ggH110;1	
 KEY: TH1D	bbH110;1	
 KEY: TH1D	ggH110_CMS_eff_t_mssmHigh_tt_13TeVDown;1	
.... 
 KEY: TH1D	ggH120;1	
 KEY: TH1D	bbH120;1	
 KEY: TH1D	ggH120_CMS_eff_t_mssmHigh_tt_13TeVDown;1	
.....
```
In the `htt_tt_9_13TeV.txt` datacard you can find the following lines: 

```
shapes * htt_tt_9_13TeV htt_input.root htt_tt_9_13TeV/$PROCESS htt_tt_9_13TeV/$PROCESS_$SYSTEMATIC
shapes bbH htt_tt_9_13TeV htt_input.root htt_tt_9_13TeV/bbH$MASS htt_tt_9_13TeV/bbH$MASS_$SYSTEMATIC
shapes ggH htt_tt_9_13TeV htt_input.root htt_tt_9_13TeV/ggH$MASS htt_tt_9_13TeV/ggH$MASS_$SYSTEMATIC

```
defining the mapping for all mass points at the same time. One can use this datacard to estimate 95%CL for different mass points by assigning <mass_value> in the command below.  

```
combine -M AsymptoticLimits  -d htt_tt_9_13TeV.txt  -m <mass_value>
```

## Simultaneous fit in multiple categories

```
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/tutorials/model_building_2024/control_regions
```

To combine the datacards corresponding to various (independent) regions into a single card one can use `combineCards.py`. 

```
combineCards.py htt_tt_9_13TeV=htt_tt_9_13TeV.txt htt_tt_8_13TeV=htt_tt_8_13TeV.txt >htt_tt_SRs.txt

```
The combined card `htt_tt_SRs.txt` now has two categories: 
```
----------------------------------------------------------------------------------------------------------------------------------
bin          htt_tt_9_13TeV  htt_tt_8_13TeV
observation  3416            105545
----------------------------------------------------------------------------------------------------------------------------------

```

## Rate parameters 

```
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/tutorials/model_building_2024/control_regions
```

It is quite common to use data-drive background estimation methods. In Combine one can perform simultaneous fit of signal and control regions it allows to automatically handle the statistical uncertainty due to the number of data events in the control region, correctly handles signal contamination in the control region, allows to properly take into account systematic uncertainties affecting the backgrounds in the control regions. 

In the working directory for this section you can find the `htt_zmm_8_13TeV.txt, htt_zmm_9_13TeV.txt and htt_ttbar_1_13TeV.txt` cards, corresponding to the control regions enriched in ZLL and ttbar processes, in addition to the signal regions from the previous step. Let's combine all of the regions into one datacard. 

```
combineCards.py htt_zmm_9_13TeV=htt_zmm_9_13TeV.txt htt_zmm_8_13TeV=htt_zmm_8_13TeV.txt htt_ttbar_1_13TeV=htt_ttbar_1_13TeV.txt htt_tt_9_13TeV=htt_tt_9_13TeV.txt htt_tt_8_13TeV=htt_tt_8_13TeV.txt >htt_tt_combined.txt
```
Now the `htt_tt_combined.txt` contains signal and control regions. To allow the rate of the background processes to be corrected from the control regions we can define common rate parameters, which linearly scale the predicted rates specified in the datacard, using the syntax 

```
<rate_param_name> rateParam <category> <process> <initial_value> [min_value,max_value]
```

The following lines define `rate_TT` and `rate_ZMM` rate parameters scaling TT and ZLL processes in all regions simultaneously: 
```
rate_TT                 rateParam  *          TT         1 [0,5]
rate_TT                 rateParam  *          TTT        1 [0,5]
rate_ZMM                rateParam  *          ZLL        1 [0,2]
```

Note that by default rate parameters are freely floating (unconstrained) parameters in Combine, however it is possible to add a constrain term to the likelihood by adding the `param` modifier with the same name as rate parameter: 

```
rate_TT param <mean> <sigma>  
```
**Task:** Add `param` nuisance named `rate_TT` with `mean = 1.` and `sigma = 1`, check how the uncertainty on `rate_TT` parameter changes, what happens if you change the width of the constraint term?
 
In addition modifiers that are functions of other parameters can be included using the following syntax:

```
name rateParam bin process formula args
```
This can be useful to constrain the ratio of two processes as shown below

```
rate_A rateParam *  process_A 
ratio_BtoA param 1 1
rate_B rateParam *  process_B  @0*@1 rate_A*ratio
```

## Extra arguments 

```
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/tutorials/model_building_2024/PhysicsModels
```

In one wants to connect different models with common parameters, or just use external functions it is possible to import parameters defined within external workspaces with `extArg`: 

```nohighlight
name extArg rootfile:workspacename
```
The `extArg` syntax allows to import `RooAbsReal` object from an external workspace. This object can be another free floating parameter, or a function of other parameters. In this section we are going to import `RooSpline1D` objects which define how various Higgs production cross sections depend on the Higgs mass (`MH` parameter). 

The datacard we are going to use in this section `htt_tt_125_8TeV.txt` correspond to 8 TeV analysis and Higgs mass of 125 GeV. To excercise the `extArg` features, let's rescale the signal templates to 13 TeV cross section values. 
The 13 and 8 TeV cross sections predictions from YR4 are stored in the `$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/`
in `sm_yr4_13TeV.root` and `sm_yr4_8TeV.root` files respectively, let's inspect the contents of `sm_yr4_13TeV.root`: 

```
 TFile*		/afs/cern.ch/work/a/anigamov/rootv630/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/sm_yr4_13TeV.root	
  KEY: RooWorkspace	xs_13TeV;1	xs_13TeV
  KEY: TProcessID	ProcessID0;1	2fd49e90-f1a0-11e8-9717-b052b8bcbeef
root [2] xs_13TeV->Print()

RooWorkspace(xs_13TeV) xs_13TeV contents

variables
---------
(MH)

functions
--------
RooSpline1D::WH_13TeV[ xvar=MH ] = 1.369
RooSpline1D::WminusH_13TeV[ xvar=MH ] = 0.5313
RooSpline1D::WplusH_13TeV[ xvar=MH ] = 0.838
RooSpline1D::ZH_13TeV[ xvar=MH ] = 0.8824
RooSpline1D::bbH_13TeV[ xvar=MH ] = 0.4863
RooSpline1D::ggH_13TeV[ xvar=MH ] = 48.52
RooSpline1D::ggZH_13TeV[ xvar=MH ] = 0.1227
RooFormulaVar::qqZH_13TeV[ actualVars=(ZH_13TeV,ggZH_13TeV) formula="@0-@1" ] = 0.7597
RooSpline1D::tHW_13TeV[ xvar=MH ] = 0.01517
RooSpline1D::tHq_13TeV[ xvar=MH ] = 0.07714
RooSpline1D::ttH_13TeV[ xvar=MH ] = 0.5065
RooSpline1D::vbfH_13TeV[ xvar=MH ] = 3.779
```
The `RooSpline1D::WH_13TeV[ xvar=MH ] = 1.369` contains cross sections values interpolated between various Higgs mass points. 
We can import them into our model as shown below

```
vbfH_13TeV     extArg     $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/sm_yr4_13TeV.root:xs_13TeV
ggH_13TeV     extArg     $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/sm_yr4_13TeV.root:xs_13TeV
ZH_13TeV     extArg     $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/sm_yr4_13TeV.root:xs_13TeV
WH_13TeV     extArg     $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/sm_yr4_13TeV.root:xs_13TeV

vbfH_8TeV     extArg     $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/sm_yr4_8TeV.root:xs_8TeV
ggH_8TeV     extArg     $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/sm_yr4_8TeV.root:xs_8TeV
ZH_8TeV     extArg     $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/sm_yr4_8TeV.root:xs_8TeV
WH_8TeV     extArg     $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/sm_yr4_8TeV.root:xs_8TeV
```
Then we can define `rateParam` as functions of imported `extArg`s to rescale signal processes WH, ZH, ggH, qqH: 
```
WH_8to13TeV    rateParam  *   WH    @0/@1         WH_13TeV,WH_8TeV
ZH_8to13TeV    rateParam  *   ZH    @0/@1         ZH_13TeV,ZH_8TeV
ggH_8to13TeV    rateParam  *   ggH    @0/@1         ggH_13TeV,ggH_8TeV
vbfH_8to13TeV    rateParam  *   qqH    @0/@1         vbfH_13TeV,vbfH_8TeV
```

When running Combine methods (e.g. `combine -M MultiDimFit --algo singles`) with these parameters you have to freeze `MH`, i.e. add `--freezeParameters MH` option.  

**Advanced task:** rescale the signal templates to the cross-section corresponding to a different MH value (e.g. 120 GeV). 

## Physics Models

With physics model one can instruct Combine how to scale the signal (background) yields with parameters of the model: $M_{cp}(\vec{\mu})$ from 

$$n^\mathrm{exp}_{cb} = \mathrm{max}(0, \sum_{p} M_{cp}(\vec{\mu})N_{cp}(\nu_G, \vec{\nu}_L,\vec{\nu}_S,\vec{\nu}_{\rho})\omega_{cbp}(\vec{\nu}_S) + E_{cb}(\vec{\nu}_B) ) $$

In Combine we can use [PhysicsModel](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/python/PhysicsModel.py#L16) as a base class.

This class has several useful methods, but the most important ones are `doParametersOfInterest()` which defines parameters and functions of the model, and `getYieldScale(self, bin, process)` defines how expected events are scaled with the model parameters.

```
class PhysicsModelBase(six.with_metaclass(ABCMeta, object)):
...
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

    def getYieldScale(self, bin, process):
        "Return the name of a RooAbsReal to scale this yield by or the two special values 1 and 0 (don't scale, and set to zero)"
        return "r" if self.DC.isSignal[process] else 1
...

```

There are many models available to use for different physics cases, follow the [link](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/part2/physicsmodels/) for more information. In the following sections we will discuss how one can construct custom models. 

To use a different physics model instead of the default one, we are going to use the option -P as in

```
text2workspace.py datacard.txt -P HiggsAnalysis.CombinedLimit.PythonFile:modelName
```

### Default physics model 

The default physics model implemented in Combine defines a single POI that linearly scales all signal processes. We use this model by default when running `text2workspace.py -m <mass_value> <datacard.txt>`. 

$$
M_{cp}(\mu) = \begin{cases}
    \mu  &\mathrm{if\ } p \in \mathrm{signal} \\
    1    &\mathrm{otherwise} \end{cases}
$$

```
class PhysicsModel(PhysicsModelBase):
    """Example class with signal strength as only POI"""

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("r[1,0,20]")
        self.modelBuilder.doSet("POI", "r")
        # --- Higgs Mass as other parameter ----
        if self.options.mass != 0:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)

```
### Multi Signal model

Combine already contains a model [`HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel`](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/part2/physicsmodels/#multisignalmodel-ready-made-model-for-multiple-signal-processes) that can be used to assign different signal strengths to multiple processes in a datacard, configurable from the command line using the mapping `--PO 'map=<bin>/<process>:<parameter_name>'`. The wildcard `*` are allowed for `<bin>` and `<process>` entries. The following command assigns `r_ggH` signal strength to scale `ggH` processes in all regions (`bin`s). 

```
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose  --PO  'map=.*/ggH:r_ggH[1,-5,5]' --PO 'map=.*/qqH:r_qqH[1,-5,5]' PhysicsModel/htt_tt_125_8TeV.txt  -o ws_multiSignal.root -m 125
combine -M MultiDimFit --algo singles -d  ws_multiSignal.root -n .multiSignal.
```

### Custom models

Now let's look at the example of how one can construct a model where the two leading signal processes (qqH and ggH) are scaled with relative fraction parameter `f` and the overall rate modifier `r`. 

$$
N_{qqH}(r, f) = r f N_{qqH};\,N_{ggH}(r, f) = r (1 - f) N_{ggH}
$$

As discussed above first we have to define the parameters of the model
```
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("f[0,0,4]") 
        self.modelBuilder.doVar("r[1,0,10]")
```
then we can use the in the scaling functions for ggH and qqH processes
```
        self.modelBuilder.factory_( "expr::scale_qqH(\"(1-@0)*@1\", f,r)")
        self.modelBuilder.factory_( "expr::scale_ggH(\"@0*@1\", f,r)")
```
add the parameters to the set of POIs
```
        self.modelBuilder.doSet("POI", ",".join(["f"]))
        self.modelBuilder.doSet("POI", ",".join(["r"]))
```
The `getYieldScale(self, bin, process)` method will scale `qqH` process with the `scale_qqH` object and `ggH` with `scale_ggH`. 
```
    def getYieldScale(self, bin, process):
        if process == "qqH": return "scale_qqH"
        elif process == "ggH": return "scale_ggH"
        else: return 1

```

To use this model we should add the directory where the corresponding file is located to the `PYTHON3PATH`:
```
export PYTHON3PATH=${PYTHON3PATH}:${PWD}/models
```
And now we can finally create the workspace using this model: 
```
text2workspace.py PhysicsModels/htt_tt_125_8TeV.txt -P FractionModel:Fraction_2signals  -m 125  -o ws_fraction.root
```

One can inspect the created workspace to ensure that the model is correctly implemented

```
root -l ws_fraction.root
root [1] w->Print()
```

The created workspace is quite large, but it should have two new `RooFormulaVar` objects 

```
RooFormulaVar::scale_ggH[ actualVars=(f,r) formula="x[0]*x[1]" ] = 0
RooFormulaVar::scale_qqH[ actualVars=(f,r) formula="(1-x[0])*x[1]" ] = 1
```
which modify the normalisation of the ggH and qqH processes

```
ProcessNormalization::n_exp_binhtt_tt_1_8TeV_proc_ggH[ thetaList=(CMS_eff_t_tt_8TeV,CMS_htt_scale_met_8TeV,QCDscale_ggH1in,UEPS,lumi_8TeV,pdf_gg) asymmThetaList=() otherFactorList=(scale_ggH) ] = 0
ProcessNormalization::n_exp_binhtt_tt_1_8TeV_proc_qqH[ thetaList=(CMS_eff_t_tt_8TeV,CMS_htt_scale_met_8TeV,CMS_scale_j_8TeV,QCDscale_qqH,UEPS,lumi_8TeV,pdf_qqbar) asymmThetaList=() otherFactorList=(scale_qqH) ] = 1.4954

```

### EFT model 

One can also define analytical BSM models with Combine. In this section we will extract the CI for one of the SMEFT Wilson coefficient. Without going into details, it can be shown that the dimension-6 SMEFT operators scale with quadratic equations. In this example will consider the $c_{Hg}$ operator that among all signal processes in the datacard affects only ggH at the LO in SMEFT. 

$\sigma(c_{g}) = \sigma_{SM} (1 + A c_{g} + B c^{2}_{g})$, where A and B coefficients are numbers that can be estimated from simulation.

```
class SMEFT_chg(PhysicsModel):
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("A[39.54]")
        self.modelBuilder.doVar("B[245.32]")
	self.modelBuilder.out.var("A").setConstant(True) 
	self.modelBuilder.out.var("B").setConstant(True)
        self.modelBuilder.doVar("chg[0,-1,1]")
        self.modelBuilder.factory_( "expr::ggH_scaling_chg(\"1+@1*@0+@2*@0*@0\", chg, A, B)")
	self.modelBuilder.doSet("POI", ",".join(["chg"]))

    def getYieldScale(self, bin, process):
        if process == "ggH": return "ggH_scaling_chg"
        else: return 1

smeft_chg_tutorial = SMEFT_chg()

```

To create the workspace using this model one can run 

```
export PYTHON3PATH=${PYTHON3PATH}:${PWD}/models
text2workspace.py PhysicsModels/htt_tt_125_8TeV.txt -P EFT_simple:smeft_chg_tutorial  -m 125  -o ws_chg.root
```

Run the likelihood scan for $c_{Hg}$ parameter and make the plot:
```
combine -M MultiDimFit --algo grid -d  ws_chg.root --setParameterRanges chg=-0.2,0.2
plot1DScan.py higgsCombineTest.MultiDimFit.mH120.root --POI chg
```
