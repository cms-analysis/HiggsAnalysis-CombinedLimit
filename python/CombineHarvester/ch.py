# This python module is essentially a thin wrapper around the real module which is embedded
# in the shared library built by scram from the C++ source code. We import everything from
# this module below, then attach a few functions that could not easily be wrapped from the
# C++, but instead are re-implemented in python in terms of other wrapped functions. The
# most notable example is the AddSyst method. The C++ version relies heavily on templates
# which is not readily adaptable to python. Instead we write the functionality entirely in a
# python function, then "attach" this function to the CombineHarvester class.

try:
    from libCombineHarvesterCombineTools import *
except ImportError:
    raise ImportError(
        'Module is missing: you need to compile the C++ -> python shared library by running make')

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
    self.ForEachObj(lambda x: res.add(func(x)))
    return res


def SetFromObs(self, func):
    res = set()
    self.ForEachObs(lambda x: res.add(func(x)))
    return res


def SetFromProcs(self, func):
    res = set()
    self.ForEachProc(lambda x: res.add(func(x)))
    return res


def SetFromSysts(self, func):
    res = set()
    self.ForEachSyst(lambda x: res.add(func(x)))
    return res

# Similar to the C++ implementation. Instead of templating on
# proxy objects that call the correct Process member functions
# just pass the function names as strings.


class SystMap:

    def __init__(self, *args):
        self.methodcallers = []
        self.tmap = dict()
        for k in args:
            self.methodcallers.append(getattr(Process, k))

    def __call__(self, *args):
        assert len(args) == len(self.methodcallers) + 1
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
    self.ForEachProc(lambda x: procs.append(x))
    if self.Verbosity() >= 1:
        print name + ':' + type
    for proc in procs:
        if not valmap.Contains(proc):
            not_added_procs.append(proc)
            continue
        tuples.discard(valmap.GetTuple(proc))
        added_procs.append(proc)
        val = valmap.GetVal(proc)
        val_d = 0.
        val_u = 0.
        formula = ''
        args = ''
        is_asymm = False
        if isinstance(val, float):
            val_u = val
        elif len(val) == 2 and isinstance(val[0], float) and isinstance(val[1], float):
            is_asymm = True
            val_d = val[0]
            val_u = val[1]
        elif len(val) == 2 and isinstance(val[0], str) and isinstance(val[1], str):
            formula = val[0]
            args = val[1]
        else:
            raise RuntimeError('Systematic value not recognised')

        target.AddSystFromProc(
            proc, name, type, is_asymm, val_u, val_d, formula, args)
    if len(tuples) > 0 and self.Verbosity() >= 1:
        print '>> Map keys that were not used to create a Systematic:'
        for tup in tuples:
            print tup
    if self.Verbosity() >= 2:
        print '>> Process entries that did not get a Systematic:'
        for proc in not_added_procs:
            print proc
        print '>> Process entries that did get a Systematic:'
        for proc in added_procs:
            print proc

# Now we turn these free functions into member functions
# of the CombineHarvester class
CombineHarvester.ParseDatacard = ParseDatacard
CombineHarvester.AddObservations = AddObservations
CombineHarvester.AddProcesses = AddProcesses
CombineHarvester.SetFromAll = SetFromAll
CombineHarvester.SetFromObs = SetFromObs
CombineHarvester.SetFromProcs = SetFromProcs
CombineHarvester.SetFromSysts = SetFromSysts
CombineHarvester.AddSyst = AddSyst
