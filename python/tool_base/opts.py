# a dictionary with some pre-defined (mostly minimizer) combine options

OPTS = {
  'vanilla'          : '--minimizerStrategy 0 --minimizerTolerance 0.1 --cminOldRobustMinimize 0',
  'prefitAsimovSToy' : '-M GenerateOnly --expectSignal 1 -t -1 --saveToys --saveWorkspace --noMCbonly 1',
  'prefitAsimovBToy' : '-M GenerateOnly --expectSignal 0 -t -1 --saveToys --saveWorkspace --noMCbonly 1',
  'robust'           : '--robustFit 1 --minimizerTolerance 0.1 --minimizerAlgo Minuit2 --minimizerStrategy 0 --minimizerAlgoForMinos Minuit2 --minimizerStrategyForMinos 0 --cminPreScan --cminPreFit 1 --X-rtd FITTER_DYN_STEP --cminFallbackAlgo "Minuit2,0:0.1" --cminFallbackAlgo "Minuit2,Minimize,0:0.1" --cminOldRobustMinimize 0',
  'robustL'          : '--robustFit 1 --minimizerTolerance 0.1 --minimizerAlgo Minuit2 --minimizerStrategy 0 --minimizerAlgoForMinos Minuit2 --minimizerStrategyForMinos 0 --cminPreScan --cminPreFit 1 --X-rtd FITTER_DYN_STEP --cminFallbackAlgo "Minuit2,0:0.1" --cminFallbackAlgo "Minuit2,Minimize,0:0.1" --cminOldRobustMinimize 0 --minimizerToleranceForMinos 0.001',
  'robustLNoScan'    : '--robustFit 1 --minimizerTolerance 0.1 --minimizerAlgo Minuit2 --minimizerStrategy 0 --minimizerAlgoForMinos Minuit2 --minimizerStrategyForMinos 0 --cminPreFit 1 --X-rtd FITTER_DYN_STEP --cminFallbackAlgo "Minuit2,0:0.1" --cminFallbackAlgo "Minuit2,Minimize,0:0.1" --cminOldRobustMinimize 0 --minimizerToleranceForMinos 0.001',
  'robustNew'        : '--robustFit 1 --minimizerTolerance 0.1 --minimizerAlgo Minuit2 --minimizerStrategy 0 --minimizerAlgoForMinos Minuit2 --minimizerStrategyForMinos 0 --cminPreScan --cminPreFit 1 --cminFallbackAlgo "Minuit2,0:0.1" --cminFallbackAlgo "Minuit2,Minimize,0:0.1" --cminOldRobustMinimize 0 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --minimizerToleranceForMinos 0.1',
  'MLHesse'          : '--minimizerTolerance 0.1 --minimizerAlgo Minuit2 --minimizerStrategy 0 --cminFallbackAlgo "Minuit2,0:0.1" --cminFallbackAlgo "Minuit2,Minimize,0:0.1" --cminOldRobustMinimize 0 --out ./ --minos none --skipBOnlyFit --noMCbonly 1 --cminPreScan'
}
