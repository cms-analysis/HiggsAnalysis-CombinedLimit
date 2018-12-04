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
```
root -ql make_MVA_8bin_ws.C
combine -M AsymptoticLimits RPV_550_TestCard.txt
combine -M Significance RPV_550_TestCard.txt -t -1 --expectSignal=1
combine -M FitDiagnostics RPV_550_TestCard.txt --plots --saveShapes --saveNormalizations
python test/diffNuisances.py --all --abs fitDiagnostics.root
root -ql fit_report_ESM.C
```

