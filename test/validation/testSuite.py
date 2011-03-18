#!/usr/bin/env python
import sys
from TestClasses import *

suite = []

for M in 'ProfileLikelihood', 'BayesianSimple', 'MarkovChainMC':
    suite += [ (M, 'fast', MultiDatacardTest(M+"_Counting", datacardGlob("simple-counting/counting-B5p5-Obs6-S*.txt"),M,"")) ]
    others = 'full' if M != 'ProfileLikelihood' else 'fast'
    suite += [ (M, others, MultiDatacardTest(M+"_Counting_Obs1" , datacardGlob("simple-counting/counting-B5p5-Obs1-S*.txt"),M,""))]
    suite += [ (M, others, MultiDatacardTest(M+"_Counting_Obs11", datacardGlob("simple-counting/counting-B5p5-Obs1-S*.txt"),M,""))]

suite += [ ('HybridNew', 'fast', MultiOptionTest("HybridNew_TestStats_Counting_Obs6", "simple-counting/counting-B5p5-Obs6-Syst30B.txt", "HybridNew",
                "--singlePoint 5 --onlyTestStat",
                {'Atlas':'--testStat=Atlas', 'LEP':'--testStat=LEP', 'TEV':'--testStat=TEV', 'Profile':'--testStat=Profile'})) ]
suite += [ ('HybridNew', 'fast', MultiOptionTest("HybridNew_TestStats_HWW", "hww4ch-1fb-B-mH140.txt", "HybridNew",
                "--singlePoint 1 --onlyTestStat",
                {'Atlas':'--testStat=Atlas', 'LEP':'--testStat=LEP', 'TEV':'--testStat=TEV'})) ]
suite += [ ('HybridNew', 'full', MultiOptionTest("HybridNew_TestStats_Counting_Obs1", "simple-counting/counting-B5p5-Obs1-Syst30B.txt", "HybridNew",
                "--singlePoint 5 --onlyTestStat",
                {'Atlas':'--testStat=Atlas', 'LEP':'--testStat=LEP', 'TEV':'--testStat=TEV', 'Profile':'--testStat=Profile'})) ]
suite += [ ('HybridNew', 'full', MultiOptionTest("HybridNew_TestStats_Counting_Obs11_S2", "simple-counting/counting-B5p5-Obs11-Syst30B.txt", "HybridNew",
                "--singlePoint 2 --onlyTestStat",
                {'Atlas':'--testStat=Atlas', 'LEP':'--testStat=LEP', 'TEV':'--testStat=TEV', 'Profile':'--testStat=Profile'})) ]
suite += [ ('HybridNew', 'full', MultiOptionTest("HybridNew_TestStats_Counting_Obs11_S10", "simple-counting/counting-B5p5-Obs11-Syst30B.txt", "HybridNew",
                "--singlePoint 10 --onlyTestStat",
                {'Atlas':'--testStat=Atlas', 'LEP':'--testStat=LEP', 'TEV':'--testStat=TEV', 'Profile':'--testStat=Profile'})) ]

suite += [ ('HybridNew', 'fast', MultiOptionTest("HybridNew_pValues_Counting_Obs1", "simple-counting/counting-B5p5-Obs6-Syst30B.txt", "HybridNew",
                "--singlePoint 5 --fork 4 -T 200 --clsAcc=1",
                {'Atlas':'--testStat=Atlas', 'LEP':'--testStat=LEP', 'TEV':'--testStat=TEV'})) ]
for X in [ "Atlas", "LEP", "TEV" ]:
    suite += [ ('HybridNew', 'fast', MultiOptionTest("HybridNew_pValues_HWW_%s_fast"%X, "hww4ch-1fb-B-mH140.txt", "HybridNew", 
                    "--singlePoint 2 --fork 2 -T 100 --clsAcc=1", {X:"--testStat="+X})) ]
    suite += [ ('HybridNew', 'full', MultiOptionTest("HybridNew_pValues_HWW_%s_full"%X, "hww4ch-1fb-B-mH140.txt", "HybridNew", 
                    "--singlePoint 2 --fork 6 -T 250 --clsAcc=1", {X:"--testStat="+X})) ]


suite += [ ('MarkovChainMC', 'full', 
            MultiOptionTest("MarkovChainMC_Proposals_Counting", "simple-counting/counting-B5p5-Obs6-Syst30U.txt", "--tries 100 -i 20000", "MarkovChainMC",
                            { 'Uniform':'--proposal=uniform', 'Gaus':'--proposal=gaus', 'Ortho':'--proposal=ortho'} )) ]

suite += [ ('ProfileLikelihood', 'fast', MultiDatacardTest("ProfileLikelihood_HWW", datacardGlob("hww4ch-1fb-B*mH*.txt"),   'ProfileLikelihood','')) ]
suite += [ ('MarkovChainMC',     'fast', MultiDatacardTest("MarkovChainMC_HWW",     datacardGlob("hww4ch-1fb-B-mH140.txt"), 'MarkovChainMC','--tries 20'))  ]
for d in datacardGlob("hww4ch-1fb-B*mH*.txt"):
    suite += [ ('MarkovChainMC',  'full',  SimpleDatacardTest(d, 'MarkovChainMC','--tries 100')) ]

suite += [ ('ProfileLikelihood', 'fast', MultiDatacardTest("ProfileLikelihood_HWW_S0", datacardGlob("hww4ch-1fb-B*mH*.txt"),   'ProfileLikelihood','-S 0')) ]
suite += [ ('MarkovChainMC',     'fast', MultiDatacardTest("MarkovChainMC_HWW_S0",     datacardGlob("hww4ch-1fb-B-mH140.txt"), 'MarkovChainMC','-S 0 --tries 50'))  ]
suite += [ ('BayesianSimple',    'fast', MultiDatacardTest("BayesianSimple_HWW_S0",    datacardGlob("hww4ch-1fb-B-mH140.txt"), 'BayesianSimple','-S 0 --tries 50')) ]


from TestSuite import *
from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options] command directory\ncommand = list, create, runLocally, runBatch, report")
parser.add_option("-M", "--method", dest="method", help="method to test", default="*", metavar="METHOD")
parser.add_option("-t", "--test",   dest="length", help="which test: fast, full", default="fast", metavar="TEST")
parser.add_option("-q", "--queue",  dest="queue",  help="queue to run in batch on", default="8nh")
parser.add_option("-f", "--format", dest="format", help="format for print output", default="text")
(options, args) = parser.parse_args()
if len(args) <= 1:
    parser.print_usage()
    sys.exit(2)

thisSuite = TestSuite(args[1], options.method, options.length, suite)
if args[0] == "list": thisSuite.listJobs()
elif args[0] == "create": thisSuite.createJobs()
elif args[0] == "runLocally": thisSuite.runLocallySync()
elif args[0] == "runBatch": thisSuite.runBatch(options.queue)
elif args[0] == "report": thisSuite.report()
elif args[0] == "print": thisSuite.printIt(options.format)
else: RuntimeError, "Unknown command %s" % args[0]
