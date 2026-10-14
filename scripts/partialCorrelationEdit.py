#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import CombineHarvester.CombineTools.ch as ch
import ROOT
import argparse
import math

# Example usage:
# python partialCorrelationEdit.py datacard.txt -m 125 \
# --process jes_1,0.5:jes_2,0.75 \
# --postfix-uncorr _2016 \
# --output-txt new_datacard.txt --output-root new_datacard.root
#
# Notes:
#  - The assumption is that the "jes_1" and "jes_2" are the names used in some other card
#    that we want to be partially correlated with the coresponding uncertainties in this card.
#    In other words: we start from the fully correlated configuration.
#  - Should only be used with lnN/shape uncertainties and non-parametric models
#  - Currently only supports positive correlations and simple decorrelation of two parameters
#  - It is essential to verify the output of this script carefully:
#      * that without supplying any --process option an identical card is produced yielding
#        the same fit result (sanity check that the CH processing was ok)
#      * after the splitting, fitting the same card should yield ~roughly the same fit result.
#        it will not be perfect because we treat all effects as if they are Guassian, which they
#        are not.

def ScaleTo(syst, val, rename=''):
    if 'shape' in syst.type():
        syst.set_scale(syst.scale() * val)
    elif 'lnN' in syst.type():
        syst.set_value_u((syst.value_u() - 1.) * val + 1.)
        if syst.asymm():
            syst.set_value_d((syst.value_d() - 1.) * val + 1.)
    else:
        raise RuntimeError('Cannot scale a systematic of type %s' % syst.type())
    if rename != '':
        syst.set_name(rename)


def Decorrelate(cb, name, correlation, postfix_corr, postfix_uncorr):
    if correlation <= 0. or correlation >= 1.:
        raise RuntimeError('Correlation coeff X must be 0 <= X < 1')
    cb_syst = cb.cp().syst_name([name])
    print('>> The following systematics will be cloned and adjusted:')
    cb_syst.PrintSysts()
    ch.CloneSysts(cb_syst, cb, lambda x: ScaleTo(x, math.sqrt(1. - correlation * correlation), name + postfix_uncorr))
    cb_syst.ForEachSyst(lambda x: ScaleTo(x, correlation, name+postfix_corr))

ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')

cb = ch.CombineHarvester()
cb.SetFlag('workspace-uuid-recycle', True)
cb.SetFlag('workspaces-use-clone', True)
cb.SetFlag('import-parameter-err', False)
cb.SetFlag("zero-negative-bins-on-import", False)
cb.SetFlag("check-negative-bins-on-import", False)

parser = argparse.ArgumentParser()

parser.add_argument('datacard', help='datacard to edit')
parser.add_argument('-m', '--mass', help='mass value', default='120')
parser.add_argument('--process', default='', help='List of nuisances to process and desired correlation coeffs: [name],[correlation]:...')
parser.add_argument('--postfix-corr', default='', help='A postfix string that will be appended to the uncertainty which is in some other card')
parser.add_argument('--postfix-uncorr', default='_uncorr', help='A postfix string that will be appended to the uncertainty which is independent of the other card')
parser.add_argument('--output-txt', default='decorrelated_card.txt')
parser.add_argument('--output-root', default='decorrelated_card.shapes.root')
args = parser.parse_args()

cb.ParseDatacard(args.datacard, mass=args.mass)

if args.process != '':
    actions = [X.split(',') for X in args.process.split(':')]

    for name, correlation in actions:
        print('>> Setting correlation coefficient of %s to %f' % (name, float(correlation)))
        if float(correlation) == 0.:
            cb.cp().RenameSystematic(cb,name,name+args.postfix_uncorr)
        else:
            Decorrelate(cb, name, float(correlation), args.postfix_corr, args.postfix_uncorr)

print('>> Writing new card and ROOT file: %s' % ((args.output_txt, args.output_root),))
cb.WriteDatacard(args.output_txt, args.output_root)
