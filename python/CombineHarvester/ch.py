# This python module is essentially a thin wrapper around the real module which is embedded
# in the shared library built by scram from the C++ source code. We import everything from
# this module below, then attach a few functions that could not easily be wrapped from the
# C++, but instead are re-implemented in python in terms of other wrapped functions. The
# most notable example is the AddSyst method. The C++ version relies heavily on templates
# which is not readily adaptable to python. Instead we write the functionality entirely in a
# python function, then "attach" this function to the CombineHarvester class.
from __future__ import absolute_import
from __future__ import print_function
import itertools
from os import environ

# Prevent cppyy's check for the PCH
environ['CLING_STANDARD_PCH'] = 'none'
import cppyy

cppyy.load_reflection_info("libCombineHarvesterCombineTools")
AutoRebin = cppyy.gbl.ch.AutoRebin
BinByBinFactory = cppyy.gbl.ch.BinByBinFactory
CardWriter = cppyy.gbl.ch.CardWriter
CheckEmptyBins = cppyy.gbl.ch.CheckEmptyBins
CheckEmptyShapes = cppyy.gbl.ch.CheckEmptyShapes
CheckNormEff = cppyy.gbl.ch.CheckNormEff
CheckSizeOfShapeEffect = cppyy.gbl.ch.CheckSizeOfShapeEffect
CheckSmallSignals = cppyy.gbl.ch.CheckSmallSignals
CloneObs = cppyy.gbl.ch.CloneObs
CloneProcs = cppyy.gbl.ch.CloneProcs
CloneProcsAndSysts = cppyy.gbl.ch.CloneProcsAndSysts
CloneSysts = cppyy.gbl.ch.CloneSysts
CombineHarvester = cppyy.gbl.ch.CombineHarvester
MassesFromRange = cppyy.gbl.ch.MassesFromRange
Object = cppyy.gbl.ch.Object
Observation = cppyy.gbl.ch.Observation
Parameter = cppyy.gbl.ch.Parameter
ParseCombineWorkspace = cppyy.gbl.ch.ParseCombineWorkspace
PrintSystematic = cppyy.gbl.ch.PrintSystematic
Process = cppyy.gbl.ch.Process
SetStandardBinNames = cppyy.gbl.ch.SetStandardBinNames
SplitSyst = cppyy.gbl.ch.SplitSyst
Systematic = cppyy.gbl.ch.Systematic
TGraphFromTable = cppyy.gbl.ch.TGraphFromTable
ValidateCards = cppyy.gbl.ch.ValidateCards
ValidateShapeTemplates = cppyy.gbl.ch.ValidateShapeTemplates
ValidateShapeUncertaintyDirection = cppyy.gbl.ch.ValidateShapeUncertaintyDirection
ValsFromRange = cppyy.gbl.ch.ValsFromRange

Observation.__str__ = Observation.to_string
Process.__str__ = Process.to_string
Systematic.__str__ = Systematic.to_string


def patch(cls):
    def patch_method(func):
        setattr(cls, "__%s__" % func.__name__, getattr(cls, func.__name__))
        setattr(cls, func.__name__, func)

    return patch_method


@patch(CombineHarvester)
def ParseDatacard(self, card, analysis="", era="", channel="", bin_id=0, mass=""):
    return self.__ParseDatacard__(card, analysis, era, channel, bin_id, mass)


@patch(CombineHarvester)
def AddObservations(self, mass=[''], analysis=[''], era=[''], channel=[''], bin=[(0, '')]):
    return self.__AddObservations__(mass, analysis, era, channel, bin)


@patch(CombineHarvester)
def AddProcesses(self, mass=[''], analysis=[''], era=[''], channel=[''], procs=[], bin=[(0, '')], signal=False):
    return self.__AddProcesses__(mass, analysis, era, channel, procs, bin, signal)


# the following wrappers use cppyy.bind_object to work around a bug
# https://github.com/root-project/root/issues/10680

@patch(CombineHarvester)
def FilterAll(self, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Object")
        return func(obj)

    self.__FilterAll__["function<bool(void*)>"](wrapper)


@patch(CombineHarvester)
def ForEachObj(self, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Object")
        func(obj)

    self.__ForEachObj__["function<void(void*)>"](wrapper)


@patch(CombineHarvester)
def FilterProcs(self, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Process")
        return func(obj)

    self.__FilterProcs__["function<bool(void*)>"](wrapper)


@patch(CombineHarvester)
def ForEachProc(self, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Process")
        func(obj)

    self.__ForEachProc__["function<void(void*)>"](wrapper)


@patch(CombineHarvester)
def FilterObs(self, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Observation")
        return func(obj)

    self.__FilterObs__["function<bool(void*)>"](wrapper)


@patch(CombineHarvester)
def ForEachObs(self, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Observation")
        func(obj)

    self.__ForEachObs__["function<void(void*)>"](wrapper)


@patch(CombineHarvester)
def FilterSysts(self, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Systematic")
        return func(obj)

    self.__FilterSysts__["function<bool(void*)>"](wrapper)


@patch(CombineHarvester)
def ForEachSyst(self, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Systematic")
        func(obj)

    self.__ForEachSyst__["function<void(void*)>"](wrapper)


def CloneObs(src, dst, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Observation")
        func(obj)

    cppyy.gbl.ch.CloneObs["function<void(void*)>"](src, dst, wrapper)


def CloneProcs(src, dst, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Process")
        func(obj)

    cppyy.gbl.ch.CloneProcs["function<void(void*)>"](src, dst, wrapper)


def CloneSysts(src, dst, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Systematic")
        func(obj)

    cppyy.gbl.ch.CloneSysts["function<void(void*)>"](src, dst, wrapper)


def CloneProcsAndSysts(src, dst, func):
    def wrapper(obj):
        obj = cppyy.bind_object(obj, "ch::Object")
        func(obj)

    cppyy.gbl.ch.CloneProcsAndSysts["function<void(void*)>"](src, dst, wrapper)


@patch(CombineHarvester)
def SetFromAll(self, func):
    res = set()
    self.ForEachObj(lambda x: res.add(func(x)))
    return res


@patch(CombineHarvester)
def SetFromObs(self, func):
    res = set()
    self.ForEachObs(lambda x: res.add(func(x)))
    return res


@patch(CombineHarvester)
def SetFromProcs(self, func):
    res = set()
    self.ForEachProc(lambda x: res.add(func(x)))
    return res


@patch(CombineHarvester)
def SetFromSysts(self, func):
    res = set()
    self.ForEachSyst(lambda x: res.add(func(x)))
    return res


# These used to produce lists in the boost::python binding

@patch(CombineHarvester)
def bin_set(self):
    return list(self.__bin_set__())


@patch(CombineHarvester)
def bin_id_set(self):
    return list(self.__bin_id_set__())


@patch(CombineHarvester)
def process_set(self):
    return list(self.__process_set__())


@patch(CombineHarvester)
def analysis_set(self):
    return list(self.__analysis_set__())


@patch(CombineHarvester)
def era_set(self):
    return list(self.__era_set__())


@patch(CombineHarvester)
def channel_set(self):
    return list(self.__channel_set__())


@patch(CombineHarvester)
def mass_set(self):
    return list(self.__mass_set__())


@patch(CombineHarvester)
def syst_name_set(self):
    return list(self.__syst_name_set__())


@patch(CombineHarvester)
def syst_type_set(self):
    return list(self.__syst_type_set__())


class SystMap:
    """
    Similar to the C++ implementation. Instead of templating on
    proxy objects that call the correct Process member functions
    just pass the function names as strings.
    """

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


@patch(CombineHarvester)
def AddSyst(self, target, name, type, valmap):
    procs = []
    not_added_procs = []
    added_procs = []
    tuples = valmap.GetTupleSet()
    self.ForEachProc(lambda x: procs.append(x))
    if self.Verbosity() >= 1:
        print(name + ':' + type)
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
        print('>> Map keys that were not used to create a Systematic:')
        for tup in tuples:
            print(tup)
    if self.Verbosity() >= 2:
        print('>> Process entries that did not get a Systematic:')
        for proc in not_added_procs:
            print(proc)
        print('>> Process entries that did get a Systematic:')
        for proc in added_procs:
            print(proc)
