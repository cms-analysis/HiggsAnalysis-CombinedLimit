import ROOT
import argparse

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input mlfit file')
parser.add_argument('-p', '--parameter', help='parameter for listing correlations')
parser.add_argument('-m', '--max', type=int, default=10, help='max number to print')

args = parser.parse_args()

fin = ROOT.TFile(args.input.split(':')[0])
rfr = fin.Get(args.input.split(':')[1])
# rfr.Print()

all_pars = rfr.floatParsFinal()
correlations = []
for i in xrange(all_pars.getSize()):
    par = all_pars.at(i)
    if par.GetName() == args.parameter:
        print '%s: %+.3f +/- %.3f' % (par.GetName(), par.getVal(), par.getError())
        continue
    correlations.append((par.GetName(),
        rfr.correlation(args.parameter, par.GetName()),
        par.getVal(),
        par.getError()))

correlations.sort(key=lambda x: abs(x[1]), reverse=True)
# print correlations

print 'covQual: %i' % rfr.covQual()
for i in xrange(min(len(correlations), args.max)):
    print '%-50s %+.3f +/- %.3f    %+.3f' % (correlations[i][0], correlations[i][2], correlations[i][3], correlations[i][1])