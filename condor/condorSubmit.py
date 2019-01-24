#!/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_3_3/external/slc6_amd64_gcc630/bin/python
####!${SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL}/external/${SCRAM_ARCH}/bin/python

import sys, os
from os import system, environ
import optparse 
import subprocess

repo = "HiggsAnalysis/CombinedLimit"

# Parse command line arguments
parser = optparse.OptionParser("usage: %prog [options]\n")

parser.add_option ('--inPut_2016',  dest='inputRoot2016',  type='string', default = 'Keras_V1.2.5_v2', help="input root file directory: 2016")
parser.add_option ('--inPut_2017',  dest='inputRoot2017',  type='string', default = 'Keras_V3.0.1_v2', help="input root file directory: 2017")

parser.add_option ('-d',  dest='signalType', type='string', default = '',             help="List of signal model, comma separated")
parser.add_option ('-t',  dest='dataType',   type='string', default = 'data',         help="Specify if running over data or sudo data")
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
x509userproxy = $ENV(X509_USER_PROXY)

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

        outputFiles = [
            "higgsCombine%s.AsymptoticLimits.mH%s.MODEL%s.root" % (options.year, mass, st),
            "higgsCombine%s%s%s.FitDiagnostics.mH%s.MODEL%s.root" % (options.year, st, mass, mass, st),
            "higgsCombine%s%s%s_SignifExp.Significance.mH%s.MODEL%s.root" % (options.year, st, mass, mass, st),
            "higgsCombineSCAN_r_wSig.MultiDimFit.mH%s.MODEL%s.root " % (mass, st),
            "higgsCombine%s.HybridNew.mH%s.MODEL%s.root" % (options.year, mass, st),
            "MVA_%s_%s_%s_ws.root"   % (options.year, st, mass),
            "MVA_2016_%s_%s_ws.root" % (st, mass),
            "MVA_2017_%s_%s_ws.root" % (st, mass),
            "ws_%s_%s_%s.root"       % (options.year, st, mass),
            "fitDiagnostics%s%s%s.root" % (options.year, st, mass), 
            "log_%s%s%s_Asymp.txt"      % (options.year, st, mass),
            "log_%s%s%s_FitDiag.txt"    % (options.year, st, mass),
            "log_%s%s%s_Sign_sig.txt"   % (options.year, st, mass),
            "log_%s%s%s_Sign_noSig.txt" % (options.year, st, mass),
            "log_%s%s%s_multiDim.txt"   % (options.year, st, mass),
            "log_%s%s%s_HybridNew.txt"  % (options.year, st, mass),
                       ]

        transfer = "transfer_output_remaps = \""
        for f in outputFiles:
            transfer += "%s = %s/output-files/%s/%s" % (f, options.outPath, outDir, f)
            if f != outputFiles[-1]:
                transfer += "; "
        transfer += "\"\n"

        fileParts.append(transfer)
        fileParts.append("Arguments = %s %s %s %s %s %s\n" % (options.inputRoot2016, options.inputRoot2017, st, mass, options.year, options.dataType))
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
