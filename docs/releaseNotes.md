# Release notes

## CMSSW 10_2_X - v8.0.0


## CMSSW 8_1_X - v7.0.13

 - Nuisance `edit` selections for bins, processes or systematic names now require a complete string match. For example, `nuisance edit add procA binA [...]` will no longer match `procAB` and `binAB`. Note that regex selections can still be used to match multiple labels, but again are now required to match the full strings.
 - Nuisance parameters can now be frozen using attributes that have been assigned to the corresponding RooRealVars. Syntax is `--freezeWithAttributes attr1,attr2,...,attrN`.
 - **[EXPERIMENTAL]** For binned analyses using autoMCStats a faster implementation of the vertical template morphing for shape uncertainties can be enabled at runtime with the option `--X-rtd FAST_VERTICAL_MORPH`. Any results using this flag should be validated carefully against the default.
 - The possibility to use some old method names has now been fully removed. When setting the `-M` option, `FitDiagnostics`, `AsymptoticLimits` and `Significance` must be used instead of, respectively, `MaxLikelihoodFit`, `Asymptotic` and `ProfileLikelihood`.
 - For Higgs analyses: added YR4 cross sections, branching ratios and partial width uncertainties in `data/lhc-hxswg/sm/`, as used in HIG-17-031
