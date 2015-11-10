import CombineHarvester.CombineTools.plotting as plot
import ROOT as R
from array import array
import json

def Tablefrom1DGraph(ta, rootfile, ma_list, filename):
  tab=int(ta)
  if (tab==0):
    return
  print "Writing table: "+filename
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

def TablefromJson(ta, jsonfile, filename):
  tab=int(ta)
  if (tab==0):
    return
  print "Writing table: "+filename
  with open(jsonfile) as json_file:    
    js = json.load(json_file)
  x = []
  jsonnames = ["-2", "-1", "expected", "+1", "+2", "observed"]
  names = ["mass", "minus2sigma", "minus1sigma", "expected", "plus1sigma", "plus2sigma", "observed"]
  
  # Get list of masses
  maxpoints = [0 for i in range(7)]
  for key in js:
    x.append(float(key))
    maxpoints[0]+=1

  # Sort list of masses (Bubblesort-Algorithm. Not very fast when there are many points to be sorted.)
  i=0
  while (i < maxpoints[0]-1):
    if (x[i+1]<x[i]):
      savex=float(x[i+1])
      x[i+1]=float(x[i])
      x[i]=float(savex)
      if (i>0):
        i-=2
    i+=1

  # Get values for the masses
  xfory = [[0.0 for i in range(6)] for j in range(maxpoints[0])] # This is incase the json-file is incomplete (values for -2,-1,exp,+1,+2,obs are missing). If everything is fine, all xfory[i] are the same as x, as well as all maxpoints.
  y = [[0.0 for i in range(6)] for j in range(maxpoints[0])]
  for i in range(maxpoints[0]):
    for k in range(6):
      if jsonnames[k] not in js[str(x[i])]:
        continue
      y[maxpoints[k+1]][k]=js[str(x[i])][jsonnames[k]]
      xfory[maxpoints[k+1]][k]=x[i]
      maxpoints[k+1]+=1

  # Prepare writing values into table
  f = open(filename, 'w')
  for name in names:
    f.write("%-20s" % name)
  f.write("\n")
  for name in names:
    f.write("--------------------")
  f.write("\n")
  
  # Calculate massrange, because range() only works for integers
  mamin = float(x[0])
  mamax = float(x[-1])
  masslist=array('d', [0]*tab)
  fl=float(mamin)
  distance=(float(mamax)-float(mamin))/float(tab-1)
  for i in range(tab):
    masslist[i]=fl
    fl+=distance

  # Calculate values and write them
  value=array('d', [0]*7)
  pas=0
  for value[0] in masslist:
    for k in range(6):
      i=1
      value[k+1]=float(-1)
      while (value[0]-xfory[i][k]>1e-11): # It's possible that the very last point in the masslist differs very slightly from mamax.
        i+=1
        if (i==maxpoints[0]): # Incase the value of -2,-1,exp,+1,+2,obs is missing for mamax.
          pas=1
          break
      if (pas==0 and value[0]>=xfory[i-1][k]): # Incase the value of -2,-1,exp,+1,+2,obs is missing for mamin.
        value[k+1]=(((y[i][k]-y[i-1][k])/(xfory[i][k]-xfory[i-1][k]))*(value[0]-xfory[i-1][k]))+y[i-1][k]
      else:
        pas=0
    for k in range(7):
      if (value[k]==float(-1)):
        f.write("%-20s" % "-----")
      else:
        f.write("%-20s" % str(value[k]))
    f.write("\n")
  f.close()
