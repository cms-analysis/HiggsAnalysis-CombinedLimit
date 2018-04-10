{
  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/interface/");
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
}
