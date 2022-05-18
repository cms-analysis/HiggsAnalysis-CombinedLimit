#!/usr/bin/env python3

from __future__ import absolute_import
from __future__ import print_function
import CombineHarvester.CombineTools.ch as ch
import ROOT as R
import glob
import numpy as np
import os
import sys
import argparse
import json
from pprint import pprint

R.PyConfig.IgnoreCommandLineOptions = True
R.gROOT.SetBatch(R.kTRUE)

R.TH1.AddDirectory(False)

parser = argparse.ArgumentParser()

parser.add_argument('cards',
                    help='Specifies the full path to the datacards to check')
parser.add_argument('--printLevel', '-p', default=1, type=int,
                    help='Specify the level of info printing (0-3, default:1)')
parser.add_argument('--readOnly', action='store_true',
                    help='If this is enabled, skip validation and only read the output json')
parser.add_argument('--checkUncertOver', '-c', default=0.1, type=float,
                    help='Report uncertainties which have a normalisation effect larger than this fraction (default:0.1)')
parser.add_argument('--reportSigUnder', '-s', default=0.001, type=float,
                    help='Report signals contributing less than this fraction of the total in a channel (default:0.001)')
parser.add_argument('--jsonFile', default='validation.json',
                    help='Path to the json file to read/write results from (default:validation.json)')
parser.add_argument('--mass', default='*',
                    help='Signal mass to use (default:*)')

args = parser.parse_args()

def print_uncertainty(js_dict, dict_key, err_type_msg):
    num_problems=0
    num_problems_peruncert = {}
    if dict_key in js_dict:
        for uncert in js_dict[dict_key]:
            probs=0;
            for mybin in js_dict[dict_key][uncert]:
                num_problems+=len(list(js_dict[dict_key][uncert][mybin].keys()))
                probs+=len(list(js_dict[dict_key][uncert][mybin].keys()))
            num_problems_peruncert[uncert]=probs
        print(">>>There were ",num_problems, "warnings of type ",err_type_msg)
        if args.printLevel > 1:
            for uncert in js_dict[dict_key]:
                print("    For uncertainty",uncert, "there were ",num_problems_peruncert[uncert]," such warnings")
                if args.printLevel > 2 :
                    for mybin in js_dict[dict_key][uncert]:
                        print("        In bin ",mybin, "the affected processes are: ",json.dumps(list(js_dict[dict_key][uncert][mybin].keys())))
    else:
        print(">>>There were no warnings of type ", err_type_msg)


def print_process(js_dict, dict_key, err_type_msg):
    num_problems=0
    num_problems_peruncert = {}
    if dict_key in js_dict:
        for mybin in js_dict[dict_key]:
            probs=0;
            num_problems+=len(js_dict[dict_key][mybin])
            num_problems_peruncert[mybin]=len(js_dict[dict_key][mybin])
        print(">>>There were ",num_problems, "warnings of type ",err_type_msg)
        if args.printLevel > 1:
            for mybin in js_dict[dict_key]:
                print("    For bin",mybin, "there were ", num_problems_peruncert[mybin]," such warnings. The affected processes are: ", json.dumps(js_dict[dict_key][mybin]))
    else:
        print(">>>There were no warnings of type", err_type_msg)

def print_process_info(js_dict, dict_key, err_type_msg):
    num_problems=0
    num_problems_peruncert = {}
    if dict_key in js_dict:
        for mybin in js_dict[dict_key]:
            probs=0;
            num_problems+=len(js_dict[dict_key][mybin])
            num_problems_peruncert[mybin]=len(js_dict[dict_key][mybin])
        print(">>>INFO: there were ",num_problems," alerts of type ",err_type_msg)
        if args.printLevel > 1:
            for mybin in js_dict[dict_key]:
                print("    For bin",mybin, "there were ", num_problems_peruncert[mybin]," such alerts. The affected processes are: ", json.dumps(js_dict[dict_key][mybin]))
    else:
        print(">>>There were no alerts of type", err_type_msg)

def print_bin(js_dict, dict_key, err_type_msg):
    num_problems=0
    num_problems_peruncert = {}
    if dict_key in js_dict:
        for mybin in js_dict[dict_key]:
            probs=0;
            num_problems+=len(js_dict[dict_key][mybin])
            num_problems_peruncert[mybin]=len(js_dict[dict_key][mybin])
        print(">>>There were ",num_problems, "warnings of type ",err_type_msg)
        if args.printLevel > 1:
            for mybin in js_dict[dict_key]:
                print("    For bin",mybin, "there were ", num_problems_peruncert[mybin]," such warnings. The affected bins of the template are: ", json.dumps(js_dict[dict_key][mybin]))
    else:
        print(">>>There were no warnings of type", err_type_msg)



cb = ch.CombineHarvester()
cb.SetFlag("check-negative-bins-on-import",0)
cb.SetFlag('workspaces-use-clone', True)

if not args.readOnly:
    cb.ParseDatacard(args.cards,"","",mass=args.mass)

    ch.ValidateCards(cb,args.jsonFile,args.checkUncertOver,args.reportSigUnder)

if args.printLevel > 0:
    print("================================")
    print("=======Validation results=======")
    print("================================")
    with open (args.jsonFile) as f:
        data = json.load(f)
        if not data:
            print(">>>There were no warnings")
        else :
            print_uncertainty(data,"uncertVarySameDirect","\'up/down templates vary the yield in the same direction\'")
            print_uncertainty(data,"uncertTemplSame","\'up/down templates are identical\'")
            print_uncertainty(data,"emptySystematicShape","\'At least one of the up/down systematic uncertainty templates is empty\'")
            print_uncertainty(data,"largeNormEff","\'Uncertainty has normalisation effect of more than %.1f%%\'"% (args.checkUncertOver*100))
            print_uncertainty(data,"smallShapeEff","\'Uncertainty probably has no genuine shape effect\'")
            print_process(data,"emptyProcessShape","\'Empty process\'")
            print_bin(data,"emptyBkgBin","\'Bins of the template empty in background\'")
            print_process_info(data,"smallSignalProc","\'Small signal process\'")


