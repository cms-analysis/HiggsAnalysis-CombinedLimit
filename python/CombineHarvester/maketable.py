import CombineHarvester.CombineTools.plotting as plot
import ROOT as R
from array import array

def Tablefrom1DGraph(ta, rootfile, ma_list, filename):
  tab=int(ta)
  mamin=ma_list[0]
  mamax=ma_list[-1]
  graph = [R.TGraph() for i in range(6)]

  # Get and sort graphs
  fin = R.TFile(rootfile, 'r')
  names = ["mass", "minus2sigma", "minus1sigma", "expected", "plus1sigma", "plus2sigma", "observed"]
  for k in range(6):
    graph[k] = plot.SortGraph(fin.Get(names[k+1]))
  fin.Close()

  # Prepare writing values into table
  f = open(filename, 'w')
  for name in names:
    f.write("%-20s" % name)
  f.write("\n")
  for name in names:
    f.write("--------------------")
  f.write("\n")

  # Calculate massrange, because range() only works for integers
  masslist=array('d', [0]*tab)
  fl=float(mamin)
  distance=(float(mamax)-float(mamin))/float(tab-1)
  for i in range(tab):
    masslist[i]=fl
    fl+=distance

  # Get values and write them
  value=array('d', [0]*7)
  for value[0] in masslist:
    for k in range(6):
      value[k+1]=graph[k].Eval(value[0])
    for k in range(7):
      f.write("%-20s" % str(value[k]))
    f.write("\n")
  f.close()
