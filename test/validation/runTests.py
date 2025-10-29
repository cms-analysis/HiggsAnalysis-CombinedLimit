#!/usr/bin/env python3

import sys
from optparse import OptionParser

import ROOT
import test_AS
import test_BS
import test_BT
import test_FC
import test_HN
import test_htt
import test_MCMC
import test_PLC
from TestClasses import *
from TestSuite import *

parser = OptionParser(usage="usage: %prog [options] command \ncommand = list, create, runLocally, runBatch, report")
parser.add_option(
    "-M",
    "--method",
    dest="method",
    help="method to test",
    default="*",
    metavar="METHOD",
)
parser.add_option(
    "-n",
    "--name",
    dest="name",
    help="name of the test directory to create; defaults to the name of the test suite",
    default=None,
    metavar="TEST",
)
parser.add_option(
    "-t",
    "--test",
    dest="suite",
    help="which test suite: fast, full",
    default="fast",
    metavar="TEST",
)
parser.add_option(
    "-s",
    "--select",
    dest="select",
    help="select these tests (regexp)",
    default=None,
    metavar="PATTERN",
)
parser.add_option(
    "-x",
    "--exclude",
    dest="exclude",
    help="exclude these tests. applied after select. ",
    default=None,
    metavar="PATTERN",
)
parser.add_option(
    "-r",
    "--reference",
    dest="reference",
    help="collate report against this one",
    default=None,
)
parser.add_option("-f", "--format", dest="format", help="format for print output", default="text")
parser.add_option(
    "-j",
    "--threads",
    dest="threads",
    help="run in parallel on N threads",
    type="int",
    default=0,
)
parser.add_option("-q", "--queue", dest="queue", help="queue to run in batch on", default="8nh")
parser.add_option(
    "-1",
    "--nofork",
    dest="nofork",
    default=False,
    action="store_true",
    help="force running on a single CPU",
)
(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    sys.exit(2)

sys.argv.append("-b-")

ROOT.gROOT.SetBatch(1)


suite = []


suite += test_AS.suite  # asymptotic CLs
suite += test_PLC.suite  # profile likelihood
suite += test_BS.suite  # bayes simple
suite += test_BT.suite  # bayes toymc
suite += test_MCMC.suite  # bayes mcmc
suite += test_HN.suite  # hybrid new
suite += test_FC.suite  # feldman cousins
suite += test_htt.suite  # all htt tests


dir = options.name if options.name else options.suite
thisSuite = TestSuite(dir, options, suite)
for cmd in args:
    if cmd == "list":
        thisSuite.listJobs()
    elif cmd == "create":
        print("Creating test suite in directory", dir)
        thisSuite.createJobs()
        print("Done.")
    elif cmd == "run":
        if options.threads:
            thisSuite.runLocallyASync(threads=options.threads)
        else:
            thisSuite.runLocallySync()
    elif cmd == "submit":
        thisSuite.runBatch(options.queue)
    elif cmd == "report":
        thisSuite.report()
        thisSuite.printIt(options.format, reference=options.reference)
    elif cmd == "print":
        thisSuite.printIt(options.format, reference=options.reference)
    else:
        RuntimeError, "Unknown command %s" % cmd
