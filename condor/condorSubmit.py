import sys, os
from os import system, environ
import optparse 
import subprocess
import time

def getOptionList(option, failMessage):
    l = []
    if option:
        l = option.split(',')
        return l
    else:
        print failMessage
        exit(0)

def main():
    repo = "HiggsAnalysis/CombinedLimit"
    seed = int(time.time())
    
    # Parse command line arguments
    parser = optparse.OptionParser("usage: %prog [options]\n")

    parser.add_option ('--inPut_2016',  dest='inputRoot2016',  type='string', default = 'Keras_V1.2.5_v2', help="input root file directory: 2016")
    parser.add_option ('--inPut_2017',  dest='inputRoot2017',  type='string', default = 'Keras_V3.0.1_v2', help="input root file directory: 2017")
    parser.add_option ('-d',            dest='signalType',     type='string', default = '',                help="List of signal model, comma separated")
    parser.add_option ('-t',            dest='dataType',       type='string', default = 'data',            help="Specify if running over data or sudo data")
    parser.add_option ('-m',            dest='masssets',       type='string', default = '',                help="List of mass models, comma separated")
    parser.add_option ('-y',            dest='year',           type='string', default = '2016',            help="year")
    parser.add_option ('-c',            dest='noSubmit', action='store_true', default = False,             help="Do not submit jobs.  Only create condor_submit.txt.")
    parser.add_option ('--output',      dest='outPath',        type='string', default = '.',               help="Name of directory where output of each condor job goes")

    parser.add_option ('--toy',         dest='toy',      action='store_true', default = False,             help="Submit toy jobs instead of the normal set of fits")
    parser.add_option ('-T',            dest='numToys',        type='int',    default = 1000,              help="Specify number of toys per job")
    parser.add_option ('--rMin',        dest='rMin',           type='float',  default = 0.05,              help="Specify minimum r value")
    parser.add_option ('--rMax',        dest='rMax',           type='float',  default = 1.00,              help="Specify maximum r value")
    parser.add_option ('--rStep',       dest='rStep',          type='float',  default = 0.05,              help="Specify step size")
    parser.add_option ('--jPerR',       dest='jPerR',          type='int',    default = 5,                 help="Specify jobs per r setting")

    options, args = parser.parse_args()
    signalType = getOptionList(options.signalType, "No dataset specified")
    masssets = getOptionList(options.masssets, "No mass model specified")
        
    executable = "run_fits.tcsh"
    if options.toy:
        executable = "run_toys.tcsh"

    fileParts = []
    fileParts.append("Universe   = vanilla\n")
    fileParts.append("Executable = %s\n" % executable)
    fileParts.append("Should_Transfer_Files = YES\n")
    fileParts.append("WhenToTransferOutput = ON_EXIT\n")
    fileParts.append("request_disk = 1000000\n")
    fileParts.append("request_memory = 4000\n")
    fileParts.append("x509userproxy = $ENV(X509_USER_PROXY)\n\n")
    
    for st in signalType:
        st = st.strip()
        for mass in masssets:
            mass = mass.strip()
            print st, mass, options.year
    
            # create the directory
            outDir = st+"_"+mass+"_"+options.year
            if not os.path.isdir("%s/output-files/%s" % (options.outPath, outDir)):
                os.makedirs("%s/output-files/%s" % (options.outPath, outDir))
    
            if not options.toy:
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
            else:
                nSteps = int(round((options.rMax - options.rMin)/options.rStep))
                for x in range(0, nSteps+1):           
                    r = options.rMin + float(x)*options.rStep 
                    print "    r = ", r
                
                    for y in range(options.jPerR):                        
                        print "        seed = ", seed

                        outputFiles = [
                            "MVA_2016_%s_%s_ws.root" % (st, mass),
                            "MVA_2017_%s_%s_ws.root" % (st, mass),
                            "ws_%s_%s_%s.root"       % (options.year, st, mass),
                            "higgsCombine%s.HybridNew.mH%s.MODEL%s.%s.root" % (options.year, mass, st, str(seed)),
                            "log_%s%s%s_%s_%s_HybridNew.txt" % (options.year, st, mass, str(r), str(seed)),
                        ]

                        transfer = "transfer_output_remaps = \""
                        for f in outputFiles:
                            transfer += "%s = %s/output-files/%s/%s" % (f, options.outPath, outDir, f)
                            if f != outputFiles[-1]:
                                transfer += "; "
                        transfer += "\"\n"

                        fileParts.append(transfer)
                        fileParts.append("Arguments = %s %s %s %s %s %s %s %s %s\n" % (options.inputRoot2016, options.inputRoot2017, st, mass, options.year, 
                                                                                       options.dataType, str(r), str(seed), str(options.numToys)))
                        fileParts.append("Output = %s/log-files/MyFit_%s_%s_%s_%s.stdout\n"%(options.outPath, st, mass, str(r), str(seed)))
                        fileParts.append("Error = %s/log-files/MyFit_%s_%s_%s_%s.stderr\n"%(options.outPath, st, mass, str(r), str(seed)))
                        fileParts.append("Log = %s/log-files/MyFit_%s_%s_%s_%s.log\n"%(options.outPath, st, mass, str(r), str(seed)))
                        fileParts.append("Queue\n\n")
                        seed+=1
    
    fout = open("condor_submit.txt", "w")
    fout.write(''.join(fileParts))
    fout.close()
    
    if not options.noSubmit: 
        system('mkdir -p %s/log-files' % options.outPath)
        system("echo 'condor_submit condor_submit.txt'")
        system('condor_submit condor_submit.txt')

if __name__ == "__main__":
    main()
