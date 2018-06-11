# Benchmark Higgs models as defined in (put ref to LHCXSWG document)

# the model equivalent to mu
from HiggsAnalysis.CombinedLimit.HiggsBenchmarkModels.CSquared import CSquaredHiggs
cSq = CSquaredHiggs()

# LHC HCG models 

# kappa models
from HiggsAnalysis.CombinedLimit.LHCHCGModels import KappaVKappaF,Kappas

cVcF = KappaVKappaF(floatbrinv=False)
cVcFinv = KappaVKappaF(floatbrinv=True)
c5 = Kappas(resolved=True)
c7 = Kappas(resolved=False)
c7inv = Kappas(resolved=False,addInvisible=True);

# lambda models 
from HiggsAnalysis.CombinedLimit.LHCHCGModels import Lambdas,LambdasReduced

lambdadu = LambdasReduced(model="ldu")
lambdalq = LambdasReduced(model="llq")
lambdafv = LambdasReduced(model="lfv")
lambda7  = Lambdas()


# Older (outdated) models ...

# Models probing the Fermion sector
from HiggsAnalysis.CombinedLimit.HiggsBenchmarkModels.FermionSectorModels import C5qlHiggs, C5udHiggs
c5ql = C5qlHiggs()
c5ud = C5udHiggs()

# Models to test Custodial symmetry
from HiggsAnalysis.CombinedLimit.HiggsBenchmarkModels.CustodialSymmetryModels import CwzHiggs, CzwHiggs, RzwHiggs, RwzHiggs, LambdaWZHiggs
lambdaWZ  = LambdaWZHiggs() 
cWZ       = CwzHiggs() 
cZW       = CzwHiggs() 
rZW       = RzwHiggs()
rWZ       = RwzHiggs()


# Minimal
from HiggsAnalysis.CombinedLimit.HiggsBenchmarkModels.MinimalModels import HiggsMinimal
higgsMinimal = HiggsMinimal()

# Model with full LO parametrization 
from HiggsAnalysis.CombinedLimit.LOFullParametrization import  PartialWidthsModel
partialWidths = PartialWidthsModel()


