This folder contains easily accessible numbers for (B)SM Higgs predictions from the [LHCHXSWG](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG).

Cross-sections and Branching ratios are sourced from the excel spreadsheet [here](https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWG/Higgs_XSBR_YR4_update.xlsx), which you will also find in this directory. 

For the coupling modifiers (kappa framework), under the `couplings` directory, the coefficients for the scaling functions are taken from the [this Twiki](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG2KAPPA).

A few important notes are below for using these numbers.

Announced in HN: https://hypernews.cern.ch/HyperNews/CMS/get/higgs/1738.html

- Uncertainties for QCD scale and PDF+a_s are given as relative uncertainties (in %)
- The ggH SM x-sections are calculated at N3LO, in these files, the uncertainty "Gauss" should be used for the scale uncertainty 
- The ZH cross-section include the contribution from ggZH but ggZH is also provided separately so that a subtraction of ZH-ggZH will give the normalisation of qqZH 
- There are numbers provided separately for the cross-sections of (W+)H and (W-)H in addition to the sum of the two (these are the LAST two columns inside 13TeV-WH.txt).
- Separate numbers for PDF and alpha_S are also provided:
- The tHW cross-sections given are repeated from the value at 125.0 GeV. This is only to allow bulding a spline.
- The cross-sections and branching ratios are also provided inside `RooWorkspaces` as splines for your convenience as a function of `MH`. You can use these inside your datacards using `rateParams` to normalise the different Higgs processes. An example datacard is provided in `data/tutorials/rate_params/simple_sm_datacard.txt`.
