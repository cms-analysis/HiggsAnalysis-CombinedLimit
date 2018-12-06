HiggsAnalysis-CombinedLimit
===========================

### Official documentation

[Manual to run combine](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/wiki)

### Standalone compilation in `lxplus`
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
source env_standalone.sh 
make -j 8; make # second make fixes compilation error of first
```

### Stealth Stop Group Setup `(LPC)`
```
bash
source /cvmfs/cms.cern.ch/cmsset_default.sh 
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_8_1_0
cd CMSSW_8_1_0/src
cmsenv
git clone git@github.com:StealthStop/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
scram b clean
scram b -j8
```

### Example

To make the RooFit workspace for combine
```
root -q -l make_MVA_8bin_ws.C
```

Calculate quick asymptotic limits.  Outputs a file called higgsCombineTest.AsymptoticLimits.mH120.root
which contains the expected (and observed) limits.
The second command does the same but instead names the file higgsCombineRPV.AsymptoticLimits.mH550.root
```
combine -M AsymptoticLimits RPV_550_TestCard.txt
combine -M AsymptoticLimits RPV_550_TestCard.txt -m 550 -n RPV
```

Calculate the significance using the asimov dataset and expected signal strength of 1
```
combine -M Significance RPV_550_TestCard.txt -t -1 --expectSignal=1
```

Run the fit and produce a root file called fitDiagnostics.root that contains RooPlots and RooFitResults:
```
combine -M FitDiagnostics RPV_550_TestCard.txt --plots --saveShapes --saveNormalizations
```

Produce formatted plots of fit results for each MVA bin:
```
root -q -l fit_report_ESM.C
```

Print the parameters resulting from the fit:
```
python test/diffNuisances.py --all --abs fitDiagnostics.root
```

To make a limit plot:
```
root -l -q makePlots.C+
```
