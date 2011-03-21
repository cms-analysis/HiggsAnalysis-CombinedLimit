from TestClasses import *

suite = []
M='Hybrid'

### Test the p-values 
for R in [ 'CLs', 'CLsplusb' ]:
    suite += [ (M, 'fast', MultiOptionTest("Counting_pValues_%s" % R, "simple-counting/counting-B5p5-Obs6-Syst30B.txt", M,
                            "--singlePoint 5 --fork 2 -T 200 --clsAcc=0.01 --rule=%s" % R,
                            {'LEP':'--testStat=LEP', 'TEV':'--testStat=TEV'})) ]
    suite += [ (M, 'full', MultiOptionTest("Counting_pValues_%s" % R, "simple-counting/counting-B5p5-Obs6-Syst30B.txt", M,
                            "--singlePoint 5 --fork 2 -T 500 --clsAcc=0.002 --rule=%s" % R,
                            {'LEP':'--testStat=LEP', 'TEV':'--testStat=TEV'})) ]

for X in [ "LEP", "TEV" ]:
        suite += [ (M, 'fast', SingleDatacardTest("HWW_pValues_%s"%X, "hww4ch-1fb-B-mH140.txt", M, 
                        "--singlePoint 2 --fork 2 -T 100 --clsAcc=1 --testStat=%s"%X )) ]
        suite += [ (M, 'full', SingleDatacardTest("HWW_pValues_%s"%X, "hww4ch-1fb-B-mH140.txt", M, 
                        "--singlePoint 2 --fork 6 -T 250 --clsAcc=1 --testStat=%s"%X )) ]

### Test the limits
optionsFast="--fork 2 -H ProfileLikelihood"
optionsFull="--fork 6 -H ProfileLikelihood --rAbsAcc=0.05 --rRelAcc=0.02"
suite += [ (M, 'fast', SingleDatacardTest("Counting", "simple-counting/counting-B5p5-Obs6-Syst30B.txt", M, options)) ]
suite += [ (M, 'fast', SingleDatacardTest("HWW",      "hww4ch-1fb-B-mH140.txt",                         M, options)) ]

suite += [ (M, 'full', MultiDatacardTest("Counting",  datacardGlob("simple-counting/counting-B5p5-Obs[16]*-S*[Uy].txt"), M, options)) ]
for d in datacardGlob("hww4ch-1fb-B*mH*.txt"):
    n = re.sub("hww4ch-1fb-(.*mH..0).txt","\\1",d)
    suite += [ (M, 'full',  SingleDatacardTest('HWW_'+n, d, M, options)) ]
