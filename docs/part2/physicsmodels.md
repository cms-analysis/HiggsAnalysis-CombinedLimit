# Physics Models

Combine can be run directly on the text based datacard. However, for more advanced physics models, the internal step to convert the datacard to a binary workspace can be performed by the user. To create a binary workspace starting from a `datacard.txt`, just do

```sh
text2workspace.py datacard.txt -o workspace.root
```

By default (without the `-o` option), the binary workspace will be named `datacard.root` - i.e the **.txt** suffix will be replaced by **.root**.

A full set of options for `text2workspace` can be found by using `--help`.

The default model which will be produced when running `text2workspace` is one in which all processes identified as signal are multiplied by a common multiplier **r**. This is all that is needed for simply setting limits or calculating significances.

`text2workspace` will convert the datacard into a pdf which summaries the analysis.
For example, lets take a look at the [data/tutorials/counting/simple-counting-experiment.txt](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/data/tutorials/counting/simple-counting-experiment.txt) datacard.

```nohighlight
# Simple counting experiment, with one signal and one background process
# Extremely simplified version of the 35/pb H->WW analysis for mH = 200 GeV,
# for 4th generation exclusion (EWK-10-009, arxiv:1102.5429v1)
imax 1  number of channels
jmax 1  number of backgrounds
kmax 2  number of nuisance parameters (sources of systematical uncertainties)
------------
# we have just one channel, in which we observe 0 events
bin         1
observation 0
------------
# now we list the expected events for signal and all backgrounds in that bin
# the second 'process' line must have a positive number for backgrounds, and 0 for signal
# then we list the independent sources of uncertainties, and give their effect (syst. error)
# on each process and bin
bin             1      1
process       ggh4G  Bckg
process         0      1
rate           4.76  1.47
------------
deltaS  lnN    1.20    -    20% uncertainty on signal
deltaB  lnN      -   1.50   50% uncertainty on background
```

If we run `text2workspace.py` on this datacard and take a look at the workspace (`w`) inside the `.root` file produced, we will find a number of different objects representing the signal, background and observed event rates as well as the nuisance parameters and signal strength **r**.

From these objects, the necessary pdf has been constructed (named `model_s`). For this counting experiment we will expect a simple pdf of the form

$$
p(n_{\mathrm{obs}}| r,\delta_{S},\delta_{B})\propto
\dfrac{[r\cdot n_{S}(\delta_{S})+n_{B}(\delta_{B})]^{n_{\mathrm{obs}}} }
{n_{\mathrm{obs}}!}e^{-[r\cdot n_{S}(\delta_{S})+n_{B}(\delta_{B})]}
\cdot e^{-\frac{1}{2}(\delta_{S}- \delta_{S}^{\mathrm{In}})^{2}}
\cdot e^{-\frac{1}{2}(\delta_{B}- \delta_{B}^{\mathrm{In}})^{2}}
$$

where the expected signal and background rates are expressed as functions of the nuisance parameters, $n_{S}(\delta_{S}) = 4.76(1+0.2)^{\delta_{S}}~$ and $~n_{B}(\delta_{B}) = 1.47(1+0.5)^{\delta_{B}}$.

The first term represents the usual Poisson expression for observing $n_{\mathrm{obs}}$ events while the second two are the Gaussian constraint terms for the nuisance parameters. In this case ${\delta^{\mathrm{In}}_S}={\delta^{\mathrm{In}}_B}=0$, and the widths of both Gaussians are 1.

A combination of counting experiments (or a binned shape datacard) will look like a product of pdfs of this kind. For a parametric/unbinned analyses, the pdf for each process in each channel is provided instead of the using the Poisson terms and a product is over the bin counts/events.

## Model building

For more complex models, `PhysicsModels` can be produced. To use a different physics model instead of the default one, use the option `-P` as in

```sh
text2workspace.py datacard -P HiggsAnalysis.CombinedLimit.PythonFile:modelName
```

Generic models can be implemented by writing a python class that:

-   defines the model parameters (by default it's just the signal strength modifier **`r`**)
-   defines how signal and background yields depend on the parameters (by default, signal scale linearly with **`r`**, backgrounds are constant)
-   potentially also modifies the systematics (e.g. switch off theory uncertainties on cross section when measuring the cross section itself)

In the case of SM-like Higgs searches the class should inherit from **`SMLikeHiggsModel`** (redefining **`getHiggsSignalYieldScale`**), while beyond that one can inherit from **`PhysicsModel`**. You can find some examples in [PhysicsModel.py](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/python/PhysicsModel.py).

In the 4-process model (`PhysicsModel:floatingXSHiggs`, you will see that each of the 4 dominant Higgs production modes get separate scaling parameters, **`r_ggH`**, **`r_qqH`**, **`r_ttH`** and **`r_VH`** (or **`r_ZH`** and **`r_WH`**) as defined in,

```python
def doParametersOfInterest(self):
  """Create POI and other parameters, and define the POI set."""
  # --- Signal Strength as only POI ---
  if "ggH" in self.modes: self.modelBuilder.doVar("r_ggH[1,%s,%s]" % (self.ggHRange[0], self.ggHRange[1]))
  if "qqH" in self.modes: self.modelBuilder.doVar("r_qqH[1,%s,%s]" % (self.qqHRange[0], self.qqHRange[1]))
  if "VH"  in self.modes: self.modelBuilder.doVar("r_VH[1,%s,%s]"  % (self.VHRange [0], self.VHRange [1]))
  if "WH"  in self.modes: self.modelBuilder.doVar("r_WH[1,%s,%s]"  % (self.WHRange [0], self.WHRange [1]))
  if "ZH"  in self.modes: self.modelBuilder.doVar("r_ZH[1,%s,%s]"  % (self.ZHRange [0], self.ZHRange [1]))
  if "ttH" in self.modes: self.modelBuilder.doVar("r_ttH[1,%s,%s]" % (self.ttHRange[0], self.ttHRange[1]))
  poi = ",".join(["r_"+m for m in self.modes])
  if self.pois: poi = self.pois
  ...
```

The mapping of which POI scales which process is handled via the following function,

```python
def getHiggsSignalYieldScale(self,production,decay, energy):
  if production == "ggH": return ("r_ggH" if "ggH" in self.modes else 1)
  if production == "qqH": return ("r_qqH" if "qqH" in self.modes else 1)
  if production == "ttH": return ("r_ttH" if "ttH" in self.modes else ("r_ggH" if self.ttHasggH else 1))
  if production in [ "WH", "ZH", "VH" ]: return ("r_VH" if "VH" in self.modes else 1)
  raise RuntimeError, "Unknown production mode '%s'" % production
```

You should note that `text2workspace` will look for the python module in `PYTHONPATH`. If you want to keep your model local, you'll need to add the location of the python file to `PYTHONPATH`.

A number of models used in the LHC Higgs combination paper can be found in [LHCHCGModels.py](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/python/LHCHCGModels.py). These can be easily accessed by providing for example `-P HiggsAnalysis.CombinedLimit.HiggsCouplings:c7` and others defined un [HiggsCouplings.py](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/python/HiggsCouplings.py).

Below are some (more generic) example models which also exist in gitHub.

### MultiSignalModel ready made model for multiple signal processes

Combine already contains a model **`HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel`** that can be used to assign different signal strengths to multiple processes in a datacard, configurable from the command line.

The model is configured passing to text2workspace one or more mappings in the form **`--PO 'map=bin/process:parameter'`**

-   **`bin`** and **`process`** can be arbitrary regular expressions matching the bin names and process names in the datacard
    Note that mappings are applied both to signals and to background processes; if a line matches multiple mappings, precedence is given to the last one in the order they are in the command line.
    it is suggested to put quotes around the argument of **`--PO`** so that the shell does not try to expand any **`*`** signs in the patterns.
-   **`parameter`** is the POI to use to scale that process (`name[starting_value,min,max]` the first time a parameter is defined, then just `name` if used more than once)
    Special values are **`1`** and **`0==; ==0`** means to drop the process completely from the card, while **`1`** means to keep the yield as is in the card with no scaling (as normally done for backgrounds); **`1`** is the default that is applied to processes that have no mappings, so it's normally not needed, but it may be used either to make the thing explicit, or to override a previous more generic match on the same command line (e.g. `--PO 'map=.*/ggH:r[1,0,5]' --PO 'map=bin37/ggH:1'` would treat ggH as signal in general, but count it as background in the channel `bin37`)

Passing the additional option **`--PO verbose`** will set the code to verbose mode, printing out the scaling factors for each process; people are encouraged to use this option to make sure that the processes are being scaled correctly.

The MultiSignalModel will define all parameters as parameters of interest, but that can be then changed from the command line of combine, as described in the following sub-section.

Some examples, taking as reference the toy datacard [test/multiDim/toy-hgg-125.txt](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/test/multiDim/toy-hgg-125.txt):

-   Scale both `ggH` and `qqH` with the same signal strength `r` (that's what the default physics model of combine does for all signals; if they all have the same systematic uncertainties, it is also equivalent to adding up their yields and writing them as a single column in the card)

```nohighlight
  $ text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/ggH:r[1,0,10]' --PO 'map=.*/qqH:r' toy-hgg-125.txt -o toy-1d.root
  [...]
  Will create a POI  r  with factory  r[1,0,10]
  Mapping  r  to  ['.*/ggH']  patterns
  Mapping  r  to  ['.*/qqH']  patterns
  [...]
  Will scale  incl/bkg  by  1
  Will scale  incl/ggH  by  r
  Will scale  incl/qqH  by  r
  Will scale  dijet/bkg  by  1
  Will scale  dijet/ggH  by  r
  Will scale  dijet/qqH  by  r
```

-   Define two independent parameters of interest `r_ggH` and `r_qqH`

```nohighlight
  $ text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/ggH:r_ggH[1,0,10]' --PO 'map=.*/qqH:r_qqH[1,0,20]' toy-hgg-125.txt -o toy-2d.root
  [...]
  Will create a POI  r_ggH  with factory  r_ggH[1,0,10]
  Mapping  r_ggH  to  ['.*/ggH']  patterns
  Will create a POI  r_qqH  with factory  r_qqH[1,0,20]
  Mapping  r_qqH  to  ['.*/qqH']  patterns
  [...]
  Will scale  incl/bkg  by  1
  Will scale  incl/ggH  by  r_ggH
  Will scale  incl/qqH  by  r_qqH
  Will scale  dijet/bkg  by  1
  Will scale  dijet/ggH  by  r_ggH
  Will scale  dijet/qqH  by  r_qqH
```

-   Fix **`ggH`** to SM, define only **`qqH`** as parameter

```nohighlight
  $ text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/ggH:1' --PO 'map=.*/qqH:r_qqH[1,0,20]' toy-hgg-125.txt -o toy-1d-qqH.root
  [...]
  Mapping  1  to  ['.*/ggH']  patterns
  Will create a POI  r_qqH  with factory  r_qqH[1,0,20]
  Mapping  r_qqH  to  ['.*/qqH']  patterns
  [...]
  Will scale  incl/bkg  by  1
  Will scale  incl/ggH  by  1
  Will scale  incl/qqH  by  r_qqH
  Will scale  dijet/bkg  by  1
  Will scale  dijet/ggH  by  1
  Will scale  dijet/qqH  by  r_qqH
```

-   Drop **`ggH`** , and define only **`qqH`** as parameter

```nohighlight
 $ text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/ggH:0' --PO 'map=.*/qqH:r_qqH[1,0,20]' toy-hgg-125.txt -o toy-1d-qqH0-only.root
 [...]
 Mapping  0  to  ['.*/ggH']  patterns
 Will create a POI  r_qqH  with factory  r_qqH[1,0,20]
 Mapping  r_qqH  to  ['.*/qqH']  patterns
 [...]
 Will scale  incl/bkg  by  1
 Will scale  incl/ggH  by  0
 Will scale  incl/qqH  by  r_qqH
 Will scale  dijet/bkg  by  1
 Will scale  dijet/ggH  by  0
 Will scale  dijet/qqH  by  r_qqH
```

### Two Hypothesis testing

The `PhysicsModel` that encodes the signal model above is the [twoHypothesisHiggs](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/python/HiggsJPC.py), which assumes that there will exist signal processes with suffix **_ALT** in the datacard. An example of such a datacard can be found under [data/benchmarks/simple-counting/twoSignals-3bin-bigBSyst.txt](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/data/benchmarks/simple-counting/twoSignals-3bin-bigBSyst.txt)

```nohighlight
 $ text2workspace.py twoSignals-3bin-bigBSyst.txt -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -m 125.7 --PO verbose -o jcp_hww.root

 MH (not there before) will be assumed to be 125.7
 Process  S  will get norm  not_x
 Process  S_ALT  will get norm  x
 Process  S  will get norm  not_x
 Process  S_ALT  will get norm  x
 Process  S  will get norm  not_x
 Process  S_ALT  will get norm  x
```

The two processes (S and S_ALT) will get different scaling parameters. The LEP-style likelihood for hypothesis testing can now be performed by setting **x** or **not_x** to 1 and 0 and comparing two likelihood evaluations.

### Signal-background interference

Since there are no such things as negative probability distribution functions, the recommended way to implement this is to start from the expression for the individual amplitudes $A$ and the parameter of interest $k$,

$$
\mathrm{Yield} = |k * A_{s} + A_{b}|^2
= k^2 * |A_{s}|^2 + k * 2 \Re(A_{s}^* A_{b}) + |A_{b}|^2
= \mu * S + \sqrt{\mu} * I + B
$$

where

$\mu = k^2, ~S = |A_{s}|^2,~B = |A_b|^2$ and $S+B+I = |A_s + A_b|^2$.

With some algebra you can work out that,

$\mathrm{Yield} = \sqrt{\mu} * \left[S+B+I\right] + (\mu-\sqrt{\mu}) * \left[S\right] + (1-\sqrt{\mu}) * \left[B\right]$

where square brackets represent the input (histograms as `TH1` or `RooDataHists`) that one needs to provide.

An example of this scheme is implemented in a [HiggsWidth](https://gitlab.cern.ch/cms-hcg/cadi/hig-17-012/-/blob/master/2l2nu/HiggsWidth.py) and is completely general, since all of the three components above are strictly positive. In this example, the POI is `CMS_zz4l_mu` and the equations for the three components are scaled (separately for the **qqH** and **ggH** processes) as,

```python
 self.modelBuilder.factory_( "expr::ggH_s_func(\"@0-sqrt(@0)\", CMS_zz4l_mu)")
 self.modelBuilder.factory_(  "expr::ggH_b_func(\"1-sqrt(@0)\", CMS_zz4l_mu)")
 self.modelBuilder.factory_(  "expr::ggH_sbi_func(\"sqrt(@0)\", CMS_zz4l_mu)")

 self.modelBuilder.factory_( "expr::qqH_s_func(\"@0-sqrt(@0)\", CMS_zz4l_mu)")
 self.modelBuilder.factory_(  "expr::qqH_b_func(\"1-sqrt(@0)\", CMS_zz4l_mu)")
 self.modelBuilder.factory_(  "expr::qqH_sbi_func(\"sqrt(@0)\", CMS_zz4l_mu)")
```

### Multi-process interference

The above formulation can be extended to multiple parameters of interest
(POIs). See
[AnalyticAnomalousCoupling](https://github.com/amassiro/AnalyticAnomalousCoupling)
for an example. However, the computational performance scales quadratically
with the number of POIs, and can get extremely expensive for 10 or more, as may
be encountered often with EFT analyses. To alleviate this issue, an accelerated
interference modeling technique is implemented for template-based analyses via
the `interferenceModel` physics model. In this model, each bin yield $y$ is parameterized

$$
y(\theta) = y_0 (\theta^\top M \theta)
$$

as a function of the POI vector $\theta$, a nominal template $y_0$, and a scaling matrix $M$.
To see how this parameterization relates to that of the previous section, we can define:

$$
y_0 = A_b^2, \qquad
M = \frac{1}{A_b^2} \begin{bmatrix}
 |A_s|^2 & \Re(A_s^* A_b) \\
 \Re(A_s A_b^*) & |A_b|^2
 \end{bmatrix}, \qquad \theta = \begin{bmatrix}
 \sqrt{\mu} \\
 1
 \end{bmatrix}
$$

which leads to the same parameterization. At present, this technique only works with
`CMSHistFunc`-based workspaces, as these are the most common workspace types
encountered and the default when using
[autoMCStats](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part2/bin-wise-stats).
To use this model, for each bin find $y_0$ and put it into the datacard as a signal process, then find $M$ and
save the lower triangular component as an array in a `scaling.json` file with a
syntax as follows:

```json
[
  {
    "channel": "my_channel",
    "process": "my_nominal_process",
    "parameters": ["sqrt_mu[1,0,2]", "Bscaling[1]"],
    "scaling": [
      [0.5, 0.1, 1.0],
      [0.6, 0.2, 1.0],
      [0.7, 0.3, 1.0]
    ]
  }
]
```

where the parameters are declared using RooFit's [factory
syntax](https://root.cern.ch/doc/v622/classRooWorkspace.html#a0ddded1d65f5c6c4732a7a3daa8d16b0)
and each row of the `scaling` field represents the scaling information of a bin, e.g. if $y_0 = |A_b|^2$
then each row would contain three entries:

$$
|A_s|^2 / |A_b|^2,\quad \Re(A_s^* A_b)/|A_b|^2,\quad 1
$$

For several coefficients, one would enumerate as follows:
```python
scaling = []
for ibin in range(nbins):
    binscaling = []
    for icoef in range(ncoef):
        for jcoef in range(icoef + 1):
            binscaling.append(amplitude_squared_for(ibin, icoef, jcoef))
    scaling.append(binscaling)
```

Then, to construct the workspace, run

```bash
text2workspace.py card.txt -P HiggsAnalysis.CombinedLimit.InterferenceModels:interferenceModel \
    --PO verbose --PO scalingData=scaling.json
```
For large amounts of scaling data, you can optionally use gzipped json (`.json.gz`) or pickle (`.pkl.gz`)
files with 2D numpy arrays for the scaling coefficients instead of lists. The function `numpy.tril_indices(ncoef)`
is helpful for extracting the lower triangle of a square matrix.

You could pick any
nominal template, and adjust the scaling as appropriate. Generally it is
advisable to use a nominal template corresponding to near where you expect the
POIs to land so that the shape systematic effects are well-modeled in that
region.

It may be the case that the relative contributions of the terms are themselves
a function of the POIs. For example, in VBF di-Higgs production, BSM
modifications to the production rate can be parameterized in the "kappa"
framework via three diagrams, with scaling coefficients $\kappa_V \kappa_\lambda$,
$\kappa_V^2$, and $\kappa_{2V}$, respectively, that interfere.  In
that case, you can declare formulas with the factory syntax to represent each
amplitude as follows:

```json
[
  {
    "channel": "a_vbf_channel",
    "process": "VBFHH",
    "parameters": ["expr::a0('@0*@1', kv[1,0,2], kl[1,0,2])", "expr::a1('@0*@0', kv[1,0,2])", "k2v[1,0,2]"],
    "scaling": [
      [3.30353674666415, -8.54170982038222, 22.96464188467882, 4.2353483207128, -11.07996258835088, 5.504469544697623],
      [2.20644332142891, -7.076836641962523, 23.50989689214267, 4.053185685866683, -13.08569222837996, 7.502346155380032]
    ]
  }
]
```

However, you will need to manually specify what the POIs should be when creating the workspace using the `POIs=` physics option, e.g.

```bash
text2workspace.py card.txt -P HiggsAnalysis.CombinedLimit.InterferenceModels:interferenceModel \
  --PO scalingData=scaling.json --PO 'POIs=kl[1,0,2]:kv[1,0,2]:k2v[1,0,2]'
```
