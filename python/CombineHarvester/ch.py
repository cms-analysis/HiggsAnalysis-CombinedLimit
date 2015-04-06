from libCHCombineTools import *
import libCHCombineTools

# a regular function
def ParseDatacard(self, card, analysis="", era="", channel="", bin_id=0, mass=""):
    return self.__ParseDatacard__(card, analysis, era, channel, bin_id, mass)

def AddObservations(self, mass=[''], analysis=[''], era=[''], channel=[''], bin=[(0, '')]):
    return self.__AddObservations__(mass, analysis, era, channel, bin)

def AddProcesses(self, mass=[''], analysis=[''], era=[''], channel=[''], procs=[], bin=[(0, '')], signal=False):
    return self.__AddProcesses__(mass, analysis, era, channel, procs, bin, signal)

def SetFromAll(self, func):
    res = set()
    self.ForEachObj(lambda x : res.add(func(x)))
    return res

def SetFromObs(self, func):
    res = set()
    self.ForEachObs(lambda x : res.add(func(x)))
    return res

def SetFromProcs(self, func):
    res = set()
    self.ForEachProc(lambda x : res.add(func(x)))
    return res

def SetFromSysts(self, func):
    res = set()
    self.ForEachSyst(lambda x : res.add(func(x)))
    return res

# now we turn it into a member function
CombineHarvester.ParseDatacard = ParseDatacard
CombineHarvester.AddObservations = AddObservations
CombineHarvester.AddProcesses = AddProcesses
CombineHarvester.SetFromAll = SetFromAll
CombineHarvester.SetFromObs = SetFromObs
CombineHarvester.SetFromProcs = SetFromProcs
CombineHarvester.SetFromSysts = SetFromSysts
