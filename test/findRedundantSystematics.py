#!/usr/bin/env python3

import json
import os
import os.path
import re
from itertools import combinations
from operator import div
from optparse import OptionParser
from pprint import pprint
from sys import argv, exit, stderr, stdout
from types import *

# search for the nuisances that have identical correlation matrix:
# given two log-normal systematics X, Y for which affect exactly the same
# set of channels and processes, they can be merged with no loss of information
# provided that
# for every channel C and process P the log(kappa(X,C,P))/log(kappa(Y,C,P)) is a constant independent of C, P
# (that's because under this condition everything only depends on a linear combination of X,Y with fixed coefficients)
# This condition that is trivially satisfied if the two systematics
# affect only a single channel and process.
# You can get the matrix by parsing the datacard in python (simplest
# example is test/datacardDump.py) and then looping on DC.systs which is
# a tuple (lsyst,nofloat,pdf,pdfargs,errline); lsyst is the name,
# errline is a map that gives you kappa[channel][process].


def filterForPDFType(allSysts, type):
    filteredSysts = [x for x in allSysts if x[2] == type]
    print("Keeping PDFs of type [%s]: kept %d of %d nuisances." % (type, len(filteredSysts), len(allSysts)))
    return filteredSysts


def all_same(items):
    return all((x - items[0]) < 1e-6 for x in items)


def asymDivide(something):
    if len(something) != 2:
        raise TypeError("asymDivision requires a pair.")

    a, b = something

    if (isinstance(a, FloatType) and isinstance(b, ListType)) or (isinstance(a, ListType) and isinstance(b, FloatType)):
        return asymDivideMixed(something)
    elif isinstance(a, ListType) and isinstance(b, ListType):
        if len(a) != len(b):
            raise TypeError("For pairs of lists, they must have the same length.")
        return asymDivideLists(something)
    else:
        raise TypeError("Don't know how to divide this data structure.")


def handleZeroes(numerator):
    if numerator:
        # if this is x/0, then it is the same as 0/x
        return 0.0
    else:
        # if this is 0/0, then return nothing as there is no correlation here
        return None


def asymDivideLists(listAndList):
    (numerators, denominators) = listAndList
    pairs = list(zip(numerators, denominators))

    quotients = list()
    for pair in pairs:
        try:
            quotients.append(pair[0] / pair[1])
        except ZeroDivisionError:
            quotients.append(handleZeroes(pair[0]))

    return quotients


def asymDivideMixed(elementAndList):
    # divide x by [y,z,...] by expanding x into [x,x,...]
    # works with both [x,[y,z,...]] and [[y,z,...],x]

    (x, yz) = elementAndList
    orderKept = True

    if isinstance(x, ListType):
        (x, yz) = (yz, x)
        orderKept = False

    xx = [x] * len(yz)

    listAndList = (xx, yz) if orderKept else (yz, xx)

    return asymDivideLists(listAndList)


def lnN_redundancies(allSysts):
    systs = filterForPDFType(allSysts, "lnN")

    systsDict = dict([(s[0], s[4]) for s in systs])

    nuisNames = [s[0] for s in systs]
    channelNames = list(systs[0][4].keys())
    processNames = list(systs[0][4][channelNames[0]].keys())

    nuisPairs = combinations(nuisNames, 2)

    #    pprint(systsDict)

    kappaRatios = dict()
    correlatedPairs = list()
    for pair in nuisPairs:
        # print 'Checking', pair
        kappaRatios[pair] = list()
        for channel in channelNames:
            for process in processNames:
                kappas = [systsDict[nuis][channel][process] for nuis in pair]

                # print 'Kappas in ', channel, process
                # pprint(kappas)

                # print 'Ratios before:'
                # pprint(kappaRatios[pair])

                try:
                    # try a simple division
                    kappaRatio = kappas[0] / kappas[1]
                    kappaRatios[pair].append(kappaRatio)
                except ZeroDivisionError:
                    kappaRatios[pair].append(handleZeroes(kappas[0]))

                except TypeError:
                    try:
                        # then this is a list of 2 values (asymmetric case)
                        kappaRatio = asymDivide(kappas)
                        kappaRatios[pair].extend(kappaRatio)
                    except TypeError as e:
                        print("Could not divide " + channel + "/" + process + ": ", e)
                        pprint(pair)
                        pprint(kappas)

                # print 'Ratios after:'
                # pprint(kappaRatios[pair])

        ### The following is an early loop termination optimization
        #                if not all_same(kappaRatios[pair]):
        #                    break
        #            if not all_same(kappaRatios[pair]):
        #            if not len(set(kappaRatios[pair])) == 1:
        #                print "Rejected ", pair
        #                print set(kappaRatios[pair])
        #                rejectedPairs.append(pair)
        #                break
        temp = [x for x in kappaRatios[pair] if x != None]
        kappaRatios[pair] = list(set(temp))
        if 0.0 not in kappaRatios[pair]:
            correlatedPairs.append(pair)
    #        print 'Ratios for ', pair
    #        pprint(kappaRatios[pair])

    #    pprint(kappaRatios)

    pprint([x for x in kappaRatios.items() if x[0] in correlatedPairs])


if __name__ == "__main__":
    # import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
    argv.append("-b-")
    import ROOT

    ROOT.gROOT.SetBatch(True)
    argv.remove("-b-")

    from HiggsAnalysis.CombinedLimit.DatacardParser import *

    parser = OptionParser(usage="usage: %prog [options] datacard.txt -o output \nrun with --help to get list of options")
    addDatacardParserOptions(parser)
    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.print_usage()
        exit(1)

    options.fileName = args[0]
    if options.fileName.endswith(".gz"):
        import gzip

        file = gzip.open(options.fileName, "rb")
        options.fileName = options.fileName[:-3]
    else:
        file = open(options.fileName, "r")

    DC = parseCard(file, options)

    lnN_redundancies(DC.systs)
