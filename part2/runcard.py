from __future__ import absolute_import

from optparse import OptionParser
from sys import exit

from HiggsAnalysis.CombinedLimit.DatacardParser import *
from HiggsAnalysis.CombinedLimit.ModelTools import *
from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.ShapeTools import *

parser = OptionParser()
addDatacardParserOptions(parser)
options, args = parser.parse_args()
options.bin = True  # make a binary workspace

DC = Datacard()
MB = None

############## Setup the datacard (must be filled in) ###########################

DC.bins = ["bin1"]  # <type 'list'>
DC.obs = {"bin1": 0.0}  # <type 'dict'>
DC.processes = ["ggH", "qqWW", "ggWW", "others"]  # <type 'list'>
DC.signals = ["ggH"]  # <type 'list'>
DC.isSignal = {
    "qqWW": False,
    "ggWW": False,
    "ggH": True,
    "others": False,
}  # <type 'dict'>
DC.keyline = [
    ("bin1", "ggH", True),
    ("bin1", "qqWW", False),
    ("bin1", "ggWW", False),
    ("bin1", "others", False),
]  # <type 'list'>
DC.exp = {"bin1": {"qqWW": 0.63, "ggWW": 0.06, "ggH": 1.47, "others": 0.22}}  # <type 'dict'>
DC.systs = [
    (
        "lumi",
        False,
        "lnN",
        [],
        {"bin1": {"qqWW": 0.0, "ggWW": 1.11, "ggH": 1.11, "others": 0.0}},
    ),
    (
        "xs_ggH",
        False,
        "lnN",
        [],
        {"bin1": {"qqWW": 0.0, "ggWW": 0.0, "ggH": 1.16, "others": 0.0}},
    ),
    (
        "WW_norm",
        False,
        "gmN",
        [4],
        {"bin1": {"qqWW": 0.16, "ggWW": 0.0, "ggH": 0.0, "others": 0.0}},
    ),
    (
        "xs_ggWW",
        False,
        "lnN",
        [],
        {"bin1": {"qqWW": 0.0, "ggWW": 1.5, "ggH": 0.0, "others": 0.0}},
    ),
    (
        "bg_others",
        False,
        "lnN",
        [],
        {"bin1": {"qqWW": 0.0, "ggWW": 0.0, "ggH": 0.0, "others": 1.3}},
    ),
]  # <type 'list'>
DC.shapeMap = {}  # <type 'dict'>
DC.hasShapes = False  # <type 'bool'>
DC.flatParamNuisances = {}  # <type 'dict'>
DC.rateParams = {}  # <type 'dict'>
DC.extArgs = {}  # <type 'dict'>
DC.rateParamsOrder = set([])  # <type 'set'>
DC.frozenNuisances = set([])  # <type 'set'>
DC.systematicsShapeMap = {}  # <type 'dict'>
DC.nuisanceEditLines = []  # <type 'list'>
DC.groups = {}  # <type 'dict'>
DC.discretes = []  # <type 'list'>


###### User defined options #############################################

options.out = "combine_workspace.root"  # Output workspace name
options.fileName = "./"  # Path to input ROOT files
options.verbose = "1"  # Verbosity

##########################################################################

if DC.hasShapes:
    MB = ShapeBuilder(DC, options)
else:
    MB = CountingModelBuilder(DC, options)

# Set physics models
MB.setPhysics(defaultModel)
MB.doModel()
