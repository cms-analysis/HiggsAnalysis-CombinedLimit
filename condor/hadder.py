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
    parser.add_option ('-d', dest='signalType',     type='string', default = '',                help="List of signal model, comma separated")
    parser.add_option ('-t', dest='dataType',       type='string', default = 'data',            help="Specify if running over data or sudo data")
    parser.add_option ('-m', dest='masssets',       type='string', default = '',                help="List of mass models, comma separated")
    parser.add_option ('-y', dest='year',           type='string', default = '2016',            help="year")
    parser.add_option ('-p', dest='inPath',         type='string', default='output-files',      help="Can pass in the input directory name")
    parser.add_option ('-o', action='store_true',                                               help="Overwrite output directory")

    options, args = parser.parse_args()    
    signalType = getOptionList(options.signalType, "No dataset specified")
    masssets = getOptionList(options.masssets, "No mass model specified")
    
    # Get input directory path
    inPath = options.inPath
    for st in signalType:
        st = st.strip()
        for mass in masssets:
            mass = mass.strip()
            print st, mass, options.year
    
            # create the directory
            outDir = st+"_"+mass+"_"+options.year
            if not os.path.isdir("%s/output-files/%s" % (options.outPath, outDir)):
                os.makedirs("%s/output-files/%s" % (options.outPath, outDir))
    
    
    # Loop over all sample options to find files to hadd
    for ds in datasets:
        #higgsCombine%s.HybridNew.mH%s.MODEL%s.%s.root" % (options.year, mass, st, str(seed))
        command = "hadd %s/%s.root %s/%s/*" % (outDir, sampleCollection, inPath, sampleCollection)
        print command
        #system(command)

if __name__ == "__main__":
    main()
