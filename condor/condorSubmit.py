#!/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_3_3/external/slc6_amd64_gcc630/bin/python
####!${SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL}/external/${SCRAM_ARCH}/bin/python

import sys, os
from os import system, environ
import optparse 
import subprocess

repo = "HiggsAnalysis/CombinedLimit"

# Parse command line arguments
parser = optparse.OptionParser("usage: %prog [options]\n")

parser.add_option ('-i',  dest='inputRoot',  type='string', default = 'Keras_V1.2.5', help="input root file directory")
parser.add_option ('-d',  dest='signalType', type='string', default = '',             help="List of signal model, comma separated")
parser.add_option ('-m',  dest='masssets',   type='string', default = '',             help="List of mass models, comma separated")
parser.add_option ('-y',  dest='year',       type='string', default = '2016',         help="year")
parser.add_option ('-c',  dest='noSubmit', action='store_true', default = False, help="Do not submit jobs.  Only create condor_submit.txt.")
parser.add_option ('--output',   dest='outPath', type='string', default = '.', help="Name of directory where output of each condor job goes")
parser.add_option ('--analyze',  dest='analyze', default = 'f', help="")

options, args = parser.parse_args()

submitFile = """Universe   = vanilla
Executable = run_fits.tcsh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
request_disk = 1000000
request_memory = 4000

"""

signalType = []
if options.signalType:
    signalType = options.signalType.split(',')
else:
    print "No dataset specified"
    exit(0)

masssets = []
if options.masssets:
    masssets = options.masssets.split(',')
else:
    print "No mass model specified"
    exit(0)

fileParts = [submitFile]
for st in signalType:
    st = st.strip()
    for mass in masssets:
        mass = mass.strip()
        print st, mass, options.year

        # create the directory
        outDir = st+"_"+mass+"_"+options.year
        if not os.path.isdir("%s/output-files/%s" % (options.outPath, outDir)):
            os.makedirs("%s/output-files/%s" % (options.outPath, outDir))

        outputFiles = ["higgsCombineTest.AsymptoticLimits.mH%s.MODEL%s.root" % (mass, st),
                       "MVA_%s_%s_%s_ws.root" % (options.year, st, mass),
                       "log_%s%s%s_Asymp.txt" % (options.year, st, mass),
                       "ws_%s_%s_%s.root"     % (options.year, st, mass),
                       ]

        #fileParts.append("transfer_output_remaps = <\"higgsCombineTest.AsymptoticLimits.mH%s.MODEL%s.root = %s/output-files/%s/higgsCombineTest.AsymptoticLimits.mH%s.MODEL%s.root;" % (mass, st, options.outPath, outDir, mass, st))
        #fileParts.append("MVA_%s_%s_%s_ws.root = %s/output-files/%s/MVA_%s_%s_%s_ws.root;" % (options.year, st, mass, options.outPath, outDir, options.year, st, mass))
        #fileParts.append("log_%s%s%s_Asymp.txt = %s/output-files/%s/log_%s%s%s_Asymp.txt;" % (options.year, st, mass, options.outPath, outDir, options.year, st, mass))
        #fileParts.append("ws_%s_%s_%s.root = %s/output-files/%s/ws_%s_%s_%s.root\"\n"      % (options.year, st, mass, options.outPath, outDir, options.year, st, mass))        

        transfer = "transfer_output_remaps = <\""
        for f in outputFiles:
            transfer += "%s = %s/output-files/%s/%s" % (f, options.outPath, outDir, f)
            if f != outputFiles[-1]:
                transfer += "; "
        transfer += "\">\n"

        fileParts.append(transfer)
        fileParts.append("Arguments = %s %s %s %s\n" % (options.inputRoot, st, mass, options.year))
        fileParts.append("Output = %s/log-files/MyFit_%s_%s.stdout\n"%(options.outPath, st, mass))
        fileParts.append("Error = %s/log-files/MyFit_%s_%s.stderr\n"%(options.outPath, st, mass))
        fileParts.append("Log = %s/log-files/MyFit_%s_%s.log\n"%(options.outPath, st, mass))
        fileParts.append("Queue\n\n")

fout = open("condor_submit.txt", "w")
fout.write(''.join(fileParts))
fout.close()

if not options.noSubmit: 
    system('mkdir -p %s/log-files' % options.outPath)
    system("echo 'condor_submit condor_submit.txt'")
    system('condor_submit condor_submit.txt')
