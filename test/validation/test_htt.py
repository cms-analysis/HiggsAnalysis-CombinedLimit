from TestClasses import *

suite = []

##############################################
M='MultiDimFit'
##############################################
minimizer='--minimizerStrategy=0 --minimizerTolerance=0.1'
suite += [ (M, 'HTT',
  WorkspaceWithExpectedTest("Summer13_HTT_CV_CF",
    datacardGlob("htt/125/htt_*_8TeV.txt"), M,
    '--default-morphing=shape2 -P "HiggsAnalysis.CombinedLimit.HiggsCouplings:cVcF" --PO "cVRange=0:4" --PO "cFRange=0:4"',
    minimizer + ' --algo grid --points 25', 125, ['CV', 'CF'])
  ) ]
suite += [ (M, 'HTT',
  WorkspaceWithExpectedTest("Summer13_HTT_RV_RF",
    datacardGlob("htt/125/htt_mt_*_8TeV.txt"), M,
    '--default-morphing=shape2 -P "HiggsAnalysis.CombinedLimit.PhysicsModel:rVrFXSHiggs"',
    minimizer + ' --algo grid --points 25', 125, ['RV', 'RF'])
  ) ]

##############################################
M='MaxLikelihoodFit'
##############################################
minimizer=('--robustFit=1 --X-rtd FITTER_NEW_CROSSING_ALGO --minimizerAlgoForMinos=Minuit2 --minimizerToleranceForMinos=0.1 '
           '--X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --minimizerAlgo=Minuit2 --minimizerStrategy=0 --minimizerTolerance=0.1')
suite += [ (M, 'HTT',
  WorkspaceWithExpectedTest("Summer13_HTT_ML", datacardGlob("htt/125/htt_mt_*_8TeV.txt"), M,
    '--default-morphing=shape2', minimizer, 125)) ]

##############################################
M='Asymptotic'
##############################################
minimizer='--minimizerStrategy=0 --minimizerTolerance=0.1'
suite += [ (M, 'HTT',
  WorkspaceWithExpectedTest("Summer13_HTT_LIMIT", datacardGlob("htt/125/htt_mt_*_8TeV.txt"), M,
    '--default-morphing=shape2', minimizer, 125)) ]

##############################################
M='ProfileLikelihood'
##############################################
minimizer=' --minimizerTolerance=0.1'
suite += [ (M, 'HTT',
  WorkspaceWithExpectedTest("Summer13_HTT_SIG_OBS", datacardGlob("htt/125/htt_mt_*_8TeV.txt"), M,
    '--default-morphing=shape2', minimizer + ' --significance' , 125)) ]
suite += [ (M, 'HTT',
  WorkspaceWithExpectedTest("Summer13_HTT_SIG_EXP", datacardGlob("htt/125/htt_mt_*_8TeV.txt"), M,
    '--default-morphing=shape2', minimizer + ' --significance --expectSignal=1 -t -1 --toysFreq' , 125)) ]
