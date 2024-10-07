# Introduction  
This file is created to keep track of issues and tests conducted for integrating AD in Combine. 
## Setup

Using daily CMSSW builds with ROOT 6.32.06: 
```
cmssw-el8 
cmsrel CMSSW_14_2_ROOT632_X_2024-10-06-2300
cd CMSSW_14_2_ROOT632_X_2024-10-06-2300/src
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git checkout roofit_ad_ichep_2024_63206_comp
cd ../../
cmsenv
scram b -j 8
```

## Tested models

### Discrete profiling 

```
 text2workspace.py data/ci/datacard_RooMultiPdf.txt.gz -o multipdf_ws.root
 python3 scripts/fitRooFitAD.py --input multipdf_ws.root
```
Fails because only `RooAbsReal` parameters are allowed, pdf indices in `RooMultiPdf` models are `RooCategory`. 
```
  RooAbsReal* RooAbsPdf::createNLL(RooAbsData& data, const RooLinkedList& cmdArgs) =>
    runtime_error: In creation of function nll_func_wrapper wrapper: input param expected to be of type RooAbsReal.
```

### autoMCStats 

```
text2workspace.py data/ci/template-analysis_shapeInterp.txt   -o template_autoMCstats_ws.root  -m 200
python3 scripts/fitRooFitAD.py --input template_autoMCstats_ws.root
```
Fails because AD is not implemented for "RooRealSumPdf" objects used to model MC statistical (https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/main/python/ShapeTools.py#L264) uncertainties with `autoMCstats` 
```
[#0] ERROR:Minimization -- An analytical integral function for class "RooRealSumPdf" has not yet been implemented.
Traceback (most recent call last):
  File "/afs/cern.ch/work/a/anigamov/CMSSW_14_2_ROOT632_X_2024-10-06-2300/src/HiggsAnalysis/CombinedLimit/scripts/fitRooFitAD.py", line 29, in <module>
    nll = pdf.createNLL(data, Constrain=constrain, GlobalObservables=global_observables, EvalBackend="codegen")
  File "/cvmfs/cms-ib.cern.ch/sw/x86_64/nweek-02858/el8_amd64_gcc12/lcg/root/6.32.06-f59fcaa47c786e7268e714e7e477ee41/lib/ROOT/_pythonization/_roofit/_rooabspdf.py", line 116, in createNLL
    return self._createNLL["RooLinkedList const&"](args[0], _pack_cmd_args(*args[1:], **kwargs))
cppyy.gbl.std.runtime_error: Could not find "createNLL<RooLinkedList const&>" (set cppyy.set_debug() for C++ errors):
  RooAbsReal* RooAbsPdf::createNLL(RooAbsData& data, const RooLinkedList& cmdArgs) =>
    runtime_error: An analytical integral function for class "RooRealSumPdf" has not yet been implemented.
```
### Template based model without MC stats  

```
text2workspace.py data/ci/template-analysis_shapeInterp_womcstats.txt   -o template_ws.root  -m 200
python3 scripts/fitRooFitAD.py --input template_ws.root
```
The fit runs with the warning shown below, but results are reasonable

```
In module 'RooFitCore':
/cvmfs/cms-ib.cern.ch/sw/x86_64/nweek-02858/el8_amd64_gcc12/lcg/root/6.32.06-f59fcaa47c786e7268e714e7e477ee41/include/RooFit/Detail/MathFuncs.h:365:45: warning: function 'LnGamma' was not differentiated because clad failed to differentiate it and no suitable overload
      was found in namespace 'custom_derivatives'
      return pdf - weight * std::log(pdf) + TMath::LnGamma(weight + 1);
                                            ^
/cvmfs/cms-ib.cern.ch/sw/x86_64/nweek-02858/el8_amd64_gcc12/lcg/root/6.32.06-f59fcaa47c786e7268e714e7e477ee41/include/RooFit/Detail/MathFuncs.h:365:45: note: falling back to numerical differentiation for 'LnGamma' since no suitable overload was found and clad
      could not derive it; to disable this feature, compile your programs with -DCLAD_NO_NUM_DIFF
```
