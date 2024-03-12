# Physics Models for Extended Higgs Sector Searches

This page lists the physics models that can be used to perform searches for additional Higgs boson measurements at the LHC, when the SM-like Higgs boson should be accounted for in the data.

## Two Higgs Models

These models are for the case where there are just Two Higgs bosons, one of which is the SM-like Higgs boson that was discovered at the LHC. The two Higgs models are implemented in the python file [`TwoHiggsModels.py`](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/python/TwoHiggsModels.py). For each of these models, we assume that the SM-like Higgs boson has mass specified by `MH_SM`, while the additional boson has mass `MH`.

You can produce the model by including the  following option in the `text2workspace.py` command:

```sh
-P HiggsAnalysis.CombinedLimit.TwoHiggsModels:model
```

| |`model`|`--PO`|POIs|<div style="width:590px">Description</div>|
| --- | ------------ | -------------------- | ---------------------- | --------------------------------------------------------------------------------------------------------- |
| Two Higgs bosons  | `twoHiggsUnconstrained`    |   N/A       | `r`, `r_SM`                   | The SM-like Higgs boson signal strength will be `r_SM` and its mass will be assumed to be `MH_SM` (default value 125.8), while the additional Higgs boson signal strength will be scaled by `r` and assumed to have mass `MH`. |
| Singlet Mixing Model | `singletMixing` |  `--PO BSMDecays`,`--PO UseVisibleMu` | `r`,`BR_BSM` | Without any options, the SM like Higgs boson will have signal strength `r`, while the additional Higgs boson is scaled by `1-r` and `BR_BSM` will not be a POI. If the option `BSMDecays` is included, the additional boson's signal strength will be `r(1-BR_BSM)` and the SM like Higgs boson will have signal strength of `1-r`. If the option `UseVisibleMu` is included too, then instead, the additional Higgs boson will get a signal strength `r` while the SM one will have `1-r(1-BR_BSM)`. |
| Singlet Mixing Model for Exclusions | `singletMixingInvisible` |  `--PO BSMDecays` | `r`,`x`,`BR_BSM` | The SM like Higgs boson will be scaled by `r*x` while the additional boson is scaled by `r*(1-x)`. If the option `BSMDecays` is included, then `BR_BSM` is also a POI and the additional Higgs boson signal strength is scaled accounting for this BSM branching fraction `r*(1-x)*(1-BR_BSM)`.|
| Two Higgs bosons with $c_{V}$, $c_{F}$ couplings | `twoHiggsCvCf` | N/A | `CV`, `CF` | Both Higgs bosons signal strengths scale accoring to the coupling to vector bosons `CV` and fermions `CF` where the scaling is determined for each contributing production/decay vertex. In this case, the additioal boson is assumed to follow the same coupling structure (i.e the SM like Higgs boson couplings). |

## Fermiophobic Higgs Model

This model is for the case where the additional Higgs boson does not couple to fermions. The fermiophobic Higgs model is implemented in the python file [`HiggsFermiophobic.py`](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/python/HiggsFermiophobic.py). In this model, the Higgs boson mass is `MH`, and it can be promoted to a POI by including the option `--PO higgsMassRange=[low,high]`.

You can produce the model by including the  following option in the `text2workspace.py` command:

```sh
-P HiggsAnalysis.CombinedLimit.HiggsFermiophobic:model_name
```

| |`model`|`--PO`|POIs|<div style="width:590px">Description</div>|
| --- | ------------ | -------------------- | ---------------------- | --------------------------------------------------------------------------------------------------------- |
| Fermiophobic Higgs  | `fp`    |   N/A       | `r`                   | The Higgs boson signal strength will be `r` for any production/decay that involves vector boson couplints only. The branching ratios are recalculated assuming no Higgs boson couplings to fermions.|
