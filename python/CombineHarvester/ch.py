from libCHCombineTools import *
import libCHCombineTools
import itertools

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

class SystMap:
  def __init__(self, *args):
      self.methodcallers = []
      self.tmap = dict()
      for k in args:
        self.methodcallers.append(getattr(Process, k))

  def __call__(self, *args):
    print len(args)
    print len(self.methodcallers)+1
    assert len(args) == len(self.methodcallers)+1
    for element in itertools.product(*(args[:-1])):
      self.tmap[element] = args[-1]
    return self

  def GetTuple(self, process):
    return tuple(x(process) for x in self.methodcallers)

  def Contains(self, process):
    return self.GetTuple(process) in self.tmap

  def GetVal(self, process):
    return self.tmap.get(self.GetTuple(process), 0.)
 
  def GetTupleSet(self):
    return set(self.tmap.keys())


def AddSyst(self, target, name, type, valmap):
  procs = []
  not_added_procs = []
  added_procs = []
  tuples = valmap.GetTupleSet()
  self.ForEachProc(lambda x : procs.append(x))
  if self.Verbosity() >= 1:
    print name + ':' + type
  for proc in procs:
    if not valmap.Contains(proc):
        not_added_procs.append(proc)
        continue
    tuples.discard(valmap.GetTuple(proc))
    added_procs.append(proc)
    val = valmap.GetVal(proc)
    is_asymm = not isinstance(val, float)
    val_d = 0.
    val_u = 0.
    if is_asymm:
      assert len(val) == 2
      val_d = val[0]
      val_u = val[1]
    else:
      val_u = val
    target.AddSystFromProc(proc, name, type, is_asymm, val_u, val_d)
  if len(tuples) > 0 and self.Verbosity() >= 1:
      print '>> Map keys that were not used to create a Systematic:'
      for tup in tuples:
        print tup
  if self.Verbosity() >= 2:
      print '>> Process entries that did not get a Systematic:'
      for proc in not_added_procs: print proc
      print '>> Process entries that did get a Systematic:'
      for proc in added_procs: print proc

# now we turn it into a member function
CombineHarvester.ParseDatacard = ParseDatacard
CombineHarvester.AddObservations = AddObservations
CombineHarvester.AddProcesses = AddProcesses
CombineHarvester.SetFromAll = SetFromAll
CombineHarvester.SetFromObs = SetFromObs
CombineHarvester.SetFromProcs = SetFromProcs
CombineHarvester.SetFromSysts = SetFromSysts
CombineHarvester.AddSyst = AddSyst
