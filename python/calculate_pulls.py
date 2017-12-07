# module for calculating various definitions of "pull" 
# each pull definition must return value, error (even if error is set to 0)
# the Asym version should return value, err hi, err lo
import sys

def allowed_methods():
  return ["relDiffAsym","compatAsym","relDiffAsymErrs"] #only allow asymm versions for now (just symmetrise errors for symmetric version)

def relDiff(x,x0,sx0):
  # Common to plot (x-x0)/sigma_{x0} - i.e difference relative to initial uncertainty
  if (sx0 == 0): return [0,1]; 
  return [(x - x0)/sx0,0]

def compat(x,x0,sx,sx0):
  # Compatibility between two extremes, value assuming no constraint, and assuming no data 
  sxo = 1./( (1./(sx*sx)) - (1./sx0*sx0))
  xo  = (x - (1./(sx0*sx0))*x0)*sxo
  return [( x - xo )/( sx*sx + sxo*sxo )**0.5,0]

def compatAsym(x,x0,sxu,sxu0,sxd,sxd0):
  # as above but we choose the uncertainties depending on whether x > x0 or x < x0
  if x > x0 : ret = compat(x,x0,sxd,sxd0)
  else : ret =  compat(x,x0,sxu,sxu0)
  ret.append(0)
  return ret

def relDiffAsym(x,x0,sxu,sxu0,sxd,sxd0):
  if x<x0: ret = relDiff(x,x0,sxd0)
  else :  ret =  relDiff(x,x0,sxu0)
  ret.append(0)
  return ret

def relDiffAsymErrs(x,x0,sxu,sxu0,sxd,sxd0):
  pull = x - x0 
  pull = (pull/sxu0) if pull >= 0 else (pull/sxd)
  pull_hi = x+sxu - x0
  pull_hi = (pull_hi/sxu0) if pull_hi >= 0 else (pull_hi/sxd0)
  pull_hi = pull_hi - pull
  pull_lo = x-sxd -x0
  pull_lo = (pull_lo/sxu0) if pull_lo >= 0 else (pull_lo/sxd0)
  pull_lo =  pull - pull_lo
  return [pull,pull_hi,pull_lo]

def returnPull(method,x,x0,sx,sx0):
  if   method == "relDiff" : return relDiff(x,x0,sx0,sx0)
  elif method == "compat"  : return compat(x,x0,sx,sx0)
  else: sys.exit("python/calculate_pulls.py -- Error, pulls method %s not understood!"%method)

def returnPullAsym(method,x,x0,sxu,sxu0,sxd,sxd0):
  if   method == "relDiffAsym" : return relDiffAsym(x,x0,sxu,sxu0,sxd,sxd0)
  elif method == "relDiffAsymErrs" : return relDiffAsymErrs(x,x0,sxu,sxu0,sxd,sxd0)
  elif method == "compatAsym"  : return compatAsym(x,x0,sxu,sxu0,sxd,sxd0)
  else: sys.exit("python/calculate_pulls.py -- Error, pulls method %s not understood!"%method)

def returnTitle(method):
  if   method in ["relDiff","relDiffAsym","relDiffAsymErrs"] : return '#hat{#theta}-#theta_{0})/#Delta#theta'
  elif method in ["compat","compatAsym"] : return '#chi^{2} compat.' 
  else: sys.exit("python/calculate_pulls.py -- Error, pulls method %s not understood!"%method)

