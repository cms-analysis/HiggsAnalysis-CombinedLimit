# Benchmark Higgs models as defined in (put ref to LHCXSWG document)

# the model equivalent to mu

from HiggsAnalysis.CombinedLimit.HiggsBenchmarkModels.CSquared import CSquaredHiggs

# Models to test Custodial symmetry
from HiggsAnalysis.CombinedLimit.HiggsBenchmarkModels.CustodialSymmetryModels import CwzHiggs, CzwHiggs, LambdaWZHiggs, RwzHiggs, RzwHiggs

# Models probing the Fermion sector
from HiggsAnalysis.CombinedLimit.HiggsBenchmarkModels.FermionSectorModels import C5qlHiggs, C5udHiggs

# Minimal
from HiggsAnalysis.CombinedLimit.HiggsBenchmarkModels.MinimalModels import HiggsMinimal

# lambda models
# kappa models
from HiggsAnalysis.CombinedLimit.LHCHCGModels import Kappas, KappaVKappaF, Lambdas, LambdasReduced

# Model with full LO parametrization
from HiggsAnalysis.CombinedLimit.LOFullParametrization import PartialWidthsModel

cSq = CSquaredHiggs()

# LHC HCG models


cVcF = KappaVKappaF(floatbrinv=False)
cVcFinv = KappaVKappaF(floatbrinv=True)
c5 = Kappas(resolved=True)
c7 = Kappas(resolved=False)
c7inv = Kappas(resolved=False, addInvisible=True)


lambdadu = LambdasReduced(model="ldu")
lambdalq = LambdasReduced(model="llq")
lambdafv = LambdasReduced(model="lfv")
lambda7 = Lambdas()


# Older (outdated) models ...


c5ql = C5qlHiggs()
c5ud = C5udHiggs()


lambdaWZ = LambdaWZHiggs()
cWZ = CwzHiggs()
cZW = CzwHiggs()
rZW = RzwHiggs()
rWZ = RwzHiggs()


higgsMinimal = HiggsMinimal()


partialWidths = PartialWidthsModel()
