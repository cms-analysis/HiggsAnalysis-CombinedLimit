import sys, os
from os import system, environ
import subprocess

import optparse
from glob import glob
import shutil

def getOptionList(option, failMessage):
    l = []
    if option:
        l = option.split(',')
        return l
    else:
        print failMessage
        exit(0)

def main():
    # Parse command line arguments
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option ('-d', dest='signalType',     type='string', default = '',     help="List of signal model, comma separated")
    parser.add_option ('-t', dest='dataType',       type='string', default = 'data', help="Specify if running over data or sudo data")
    parser.add_option ('-m', dest='masssets',       type='string', default = '',     help="List of mass models, comma separated")
    parser.add_option ('-y', dest='year',           type='string', default = '2016', help="year")
    parser.add_option ('-p', dest='inPath',         type='string', default ='??',    help="Can pass in the input directory name")
    parser.add_option ('-f', action='store_true',                                    help="Overwrite output directory")

    options, args = parser.parse_args()    
    signalType = getOptionList(options.signalType, "No dataset specified")
    masssets = getOptionList(options.masssets, "No mass model specified")
    overwrite = '-f' if options.f else ''

    inPath = options.inPath
    for st in signalType:
        st = st.strip()
        for mass in masssets:
            mass = mass.strip()
            print "-----------------------------------------------------------"
            print st, mass, options.year
            print "-----------------------------------------------------------"

            path = options.inPath+"/output-files/"+st+"_"+mass+"_"+options.year
            infiles = "higgsCombine%s.HybridNew.mH%s.MODEL%s.*.root" % (options.year, mass, st)
            haddOutfile = "toys_%s_%s_%s_data.root" % (st, mass, options.year)
            command = "hadd %s %s/%s %s/%s" % (overwrite, path, haddOutfile, path, infiles)
            system(command)

            workSpace = "%s/ws_%s_%s_%s.root"%(path,options.year,st,mass)
            baseCommand = "combine -M HybridNew --LHCmode LHC-limits {0} -m {1} --keyword-value MODEL={2} -n {3} --readHybridResults --grid={4}/{5} --fullGrid".format(workSpace, mass, st, options.year, path, haddOutfile)
            commands = []
            commands.append(baseCommand + " --plot=%s/limit_scan.png" % path)
            commands.append("mv higgsCombine{0}.HybridNew.mH{1}.MODEL{2}.root {3}/higgsCombine{0}.HybridNew.mH{1}.MODEL{2}.observedLi.root".format(options.year,mass,st,path))

            for cl in ["0.025","0.160","0.500","0.840","0.975"]:
                commands.append(baseCommand + " --plot={0}/limit_scan_expected.quant{1}.png --expectedFromGrid={1}".format(path,cl))
                commands.append("mv higgsCombine{0}.HybridNew.mH{1}.MODEL{2}.quant{4}.root {3}/.".format(options.year,mass,st,path,cl))

            for c in commands:
                system(c)
            
if __name__ == "__main__":
    main()
