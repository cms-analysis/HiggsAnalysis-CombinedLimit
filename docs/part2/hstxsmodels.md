# Physics Models for STXS to EFT Interpretations

These physics models interpret Simplified Template Cross Section (STXS) measurements in terms of dimension-six effective field theory (EFT) coefficients. They scale Higgs production bins and decay branching ratios using the parameterisations supplied with Combine. Use signal process names containing an STXS production bin and a Higgs decay suffix (for example, `_hgg` or `_hzz`).

## Higgs Effective Lagrangian: stages 0, 1 and 1.1

The models in [`STXStoEFTModel.py`](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/python/STXStoEFTModel.py) use Higgs Effective Lagrangian (HEL) [ [JHEP 1307 (2013) 035](
https://link.springer.com/article/10.1007/JHEP07(2013)035), [JHEP 1404 (2014) 110](https://link.springer.com/article/10.1007/JHEP04(2014)110) ] parameterisations for stage 0, stage 1 and stage 1.1 STXS bins. The coefficient ranges and production and decay scaling functions are read from the [HEL input files](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/tree/main/data/lhc-hxswg/eft/HEL). Select a model with:

```sh
text2workspace.py datacard.txt -P HiggsAnalysis.CombinedLimit.STXStoEFTModel:model
```

| <div style="width:120px"></div> | <div style="width:100px"></div> `model` | <div style="width:190px">`--PO`</div> | <div style="width:220px">POIs</div> | <div style="width:460px">Description</div> |
| --- | --- | --- | --- | --- |
| Stage 0 | `Stage0toEFT` | `--PO freezeOtherParameters=0`, `--PO BRU=1`, `--PO STXSU=1`, `--PO fixProcesses=bin1,bin2`, `--PO higgsMassRange=x,y` | HEL coefficients (see below) | Scales stage 0 production bins and Higgs branching ratios. |
| Stage 1 | `Stage1toEFT` | Same as stage 0 | HEL coefficients (see below) | Scales stage 1 production bins and Higgs branching ratios. |
| Stage 1.1 | `Stage1_1toEFT` | Same as stage 0 | HEL coefficients (see below) | Scales stage 1.1 production bins and Higgs branching ratios. |
| Combined stages | `AllStagestoEFT` | Same as stage 0; also `--PO linearOnly=1`, `--PO useExtendedVBFScheme=1`, `--PO useLHCHXSWGStage1=1` | HEL coefficients (see below) | Accepts stage 0, 1 and 1.1 bins in one model; uses the most recent available stage for bins shared between stages. |

By default, the HEL models fit `cG_x05`, `cA_x04`, `cu_x01`, `cd_x01`, `cl_x01`, `cHW_x02`, `cHB_x01` and `cWWMinuscB_x02`. The suffixes encode the coefficient rescaling used by the model (for example, `_x02` represents a factor of $10^{-2}$). `cWWPluscB_x02` is fixed to zero by default but is used to construct the `cWW` and `cB` combinations. 

  * Set `--PO freezeOtherParameters=0` to include all coefficients listed in the HEL `pois.txt` input as POIs. 
  * `--PO fixProcesses=...` leaves the named production bins at their SM yield instead of applying EFT scaling.
  * `--PO BRU=1` adds partial-width theory uncertainties. 
  * `--PO STXSU=1` enables STXS-bin theory uncertainties, but the model warns that the supplied bin uncertainties need updating; use them with caution. 
  * For the combined-stage model, `linearOnly=1` keeps only terms linear in the coefficients, `useExtendedVBFScheme=1` selects the extended VBF production inputs, and `useLHCHXSWGStage1=1` selects the alternative stage 1 coefficients. 
  * `--PO higgsMassRange=x,y` floats the Higgs mass between $x$ and $y$; the parameterisations assume $m_H=125$ GeV, so the model warns when this option is used.

## SMEFT: stage 1.2

The Stage 1.2 interpretation uses the SMEFT [[Phys. Rept. 793 (2019) 1-98](https://www.sciencedirect.com/science/article/abs/pii/S0370157318303223?via%3Dihub)], and is implemented in [`STXStoSMEFTModel.py`](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/python/STXStoSMEFTModel.py). It reads Wilson-coefficient definitions from `pois.yaml` and production and decay terms from `prod.json` and `decay.json` in the selected [SMEFT parameterisation directory](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/tree/main/data/eft/STXStoSMEFT). Select it with:

```sh
text2workspace.py datacard.txt -P HiggsAnalysis.CombinedLimit.STXStoSMEFTModel:STXStoSMEFT
```

| <div style="width:120px"></div> | <div style="width:100px"></div>  `model` | <div style="width:190px">`--PO`</div> | <div style="width:220px">POIs</div> | <div style="width:460px">Description</div> |
| --- | --- | --- | --- | --- |
| Stage 1.2 SMEFT | `STXStoSMEFT` | `--PO parametrisation=name`, `--PO linear_only=1`, `--PO linquad_only=1`, `--PO expand_equations=1`, `--PO stage0=1`, `--PO eigenvalueThreshold=value`, `--PO fixProcesses=bin1,bin2`, `--PO higgsMassRange=x,y` | Wilson coefficients from the selected `pois.yaml`; see below | Scales production and decay using SMEFT interference and quadratic terms. |

The default `parametrisation` is `CMS-prelim-SMEFT-topU3l_22_05_05_AccCorr_0p01`, the parameterisation shipped with the repository. Its `pois.yaml` defines the following POIs (including the rescaling suffixes used in the workspace):

- Positive exponents: `chgXE3`, `chwXE2`, `chbXE2`, `chwbXE2`, `cbhreXE2`, `chj1XE1`, `chj3XE1`, `chuXE1`, `chdXE1`, `ctgreXE1`, `ctwreXE1`, `ctbreXE2`, `cwXE1`.
- Negative exponents: `chtXEm2`, `chtbreXEm1`, `cuhreXEm2`.
- No suffix: `chgtil`, `chbox`, `chdd`, `cbgre`, `chq1`, `chq3`, `chl3`, `cll1`, `chwtil`, `cbgim`, `cbhim`, `cgtil`, `cg`, `chwbtil`, `chbtil`, `cthre`, `chl1`, `che`, `chbq`, `cqj31`, `cqj38`, `ctgim`, `ctwim`, `cbwim`, `chtbim`, `cbwre`, `cthim`, `cqj18`, `ctu8`, `ctd8`, `cqu8`, `ctj8`, `cqd8`, `cqj11`, `ctu1`, `ctd1`, `cqu1`, `ctj1`, `cqd1`, `ctbim`, `cehre`.

A nonzero `exponent` adds an `XE` or `XEm` suffix to the coefficient name to account for its rescaling. 

   * By default the workspace also contains a separate set of linear-only parameters prefixed with `l`, while the regular POIs control the linear-plus-quadratic scaling. Use `linear_only=1` or `linquad_only=1` to build just one of these sets. 
   * `expand_equations=1` uses a Taylor expansion of the combined production-times-branching-ratio scaling rather than the default separate production, partial-width and total-width scaling functions.

   * `stage0=1` groups STXS bins by production mode before looking up their scaling. 
   * `eigenvalueThreshold=value` retains coefficients without an `eigenvalue` entry and those whose eigenvalue exceeds the threshold; the default is `-1` (no filtering). 
   * `fixProcesses=...` keeps the listed production bins at their SM yield. 
   * As in the HEL models, `higgsMassRange=x,y` floats the mass despite the parameterisation's $m_H=125$ GeV assumption.

