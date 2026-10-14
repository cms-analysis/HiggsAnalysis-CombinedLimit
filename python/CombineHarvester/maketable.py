from __future__ import absolute_import
import CombineHarvester.CombineTools.plotting as plot
import ROOT as R
from array import array
import json
from six.moves import range

def Tablefrom1DGraph(rootfile, filename):
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

    # Get values and write them
    value=array('d', [0]*7)
    for i in range(graph[1].GetN()) :
        value[0] = graph[1].GetX()[i]
        for k in range(6):
            value[k+1]=graph[k].GetY()[i]
        for k in range(7):
            f.write("%-20s" % str(value[k]))
        f.write("\n")
    f.close()


def TablefromJson(jsonfile, filename):
    with open(jsonfile) as json_file:
        js = json.load(json_file)
    x = []
    jsonnames = ["-2", "-1", "expected", "+1", "+2", "observed"]
    names = ["mass", "minus2sigma", "minus1sigma", "expected", "plus1sigma", "plus2sigma", "observed"]

    # Get list of masses
    maxpoints = 0
    for key in js:
        x.append(float(key))
        maxpoints+=1

    # Sort list of masses (Bubblesort-Algorithm. Not very fast when there are many points to be sorted.)
    i=0
    while (i < maxpoints-1):
        if (x[i+1]<x[i]):
            savex=float(x[i+1])
            x[i+1]=float(x[i])
            x[i]=float(savex)
            if (i>0):
                i-=2
        i+=1

    # Get values for the masses
    y = [[0.0 for i in range(6)] for j in range(maxpoints)]
    for i in range(maxpoints):
        for k in range(6):
            if jsonnames[k] in js[str(x[i])]:
                y[i][k]=js[str(x[i])][jsonnames[k]]

    # Prepare writing values into table
    f = open(filename, 'w')
    for name in names:
        f.write("%-20s" % name)
    f.write("\n")
    for name in names:
        f.write("--------------------")
    f.write("\n")

    # Write table
    i=0
    for mass in x :
        for k in range(6):
            if k==0:
                f.write("%-20s" % str(mass))
            f.write("%-20s" % str(y[i][k]))
        f.write("\n")
        i+=1
    f.close()
