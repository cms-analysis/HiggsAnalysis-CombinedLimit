# module for calculating various definitions of "pull" 
# each pull definition must return value, error (even if error is set to 0)
# the Asym version should return value, err hi, err lo
import sys

def allowed_methods():
  return ["unconstPullAsym","compatAsym","relDiffAsymErrs","diffPullAsym"] #only allow asymm versions for now (just symmetrise errors for symmetric version)

def unconstPull(x,x0,sx):
  # Common to plot (x-x0)/sigma - i.e difference relative to post-fit uncertainty
  if (sx == 0): return [0,1]; 
  return [(x - x0)/sx,0]

def compat(x,x0,sx,sx0):
  # Compatibility between two extremes, value assuming no constraint, and assuming no data 
  sxo = 1./( (1./(sx*sx)) - (1./sx0*sx0))
  xo  = (x - (1./(sx0*sx0))*x0)*sxo
  return [( x - xo )/( sx*sx + sxo*sxo )**0.5,0]

def diffPull(x,x0,sx,sx0):
  # as defined in http://physics.rockefeller.edu/luc/technical_reports/cdf5776_pulls.pdf
  if abs(sx*sx - sx0*sx0) < 0.001: 
     return [0,999]
  elif sx > sx0: 
     return [0,999]
  else:
     return [( x - x0 )/( sx0*sx0 - sx*sx )**0.5,0]
    

def compatAsym(x,x0,sxu,sxu0,sxd,sxd0):
  # as above but we choose the uncertainties depending on whether x > x0 or x < x0
  if x > x0 : ret = compat(x,x0,sxd,sxd0)
  else : ret =  compat(x,x0,sxu,sxu0)
  ret.append(0)
  return ret

def diffPullAsym(x,x0,sxu,sxu0,sxd,sxd0):
  # as defined in http://physics.rockefeller.edu/luc/technical_reports/cdf5776_pulls.pdf
  if x > x0 : ret = diffPull(x,x0,sxd,sxd0)
  else : ret =  diffPull(x,x0,sxu,sxu0)
  ret.append(ret[-1])
  return ret

def unconstPullAsym(x,x0,sxu,sxu0,sxd,sxd0):
  if x<x0: ret = unconstPull(x,x0,sxu)
  else :  ret =  unconstPull(x,x0,sxd)
  ret.append(0)
  return ret

def relDiffAsymErrs(x,x0,sxu,sxu0,sxd,sxd0):
  pull = x - x0 
  pull = (pull/sxu0) if pull >= 0 else (pull/sxd0)
  pull_hi = x+sxu - x0
  pull_hi = (pull_hi/sxu0) if pull_hi >= 0 else (pull_hi/sxd0)
  pull_hi = pull_hi - pull
  pull_lo = x-sxd -x0
  pull_lo = (pull_lo/sxu0) if pull_lo >= 0 else (pull_lo/sxd0)
  pull_lo =  pull - pull_lo
  return [pull,pull_hi,pull_lo]

def returnPull(method,x,x0,sx,sx0):
  if   method == "unconstPull" : return unconstPullDiff(x,x0,sx,sx)
  elif method == "compat"  : return compat(x,x0,sx,sx0)
  elif method == "diffPull"  : return diffPull(x,x0,sx,sx0)
  else: sys.exit("python/calculate_pulls.py -- Error, pulls method %s not understood!"%method)

def returnPullAsym(method,x,x0,sxu,sxu0,sxd,sxd0):
  if   method == "unconstPullAsym"   : return unconstPullAsym(x,x0,sxu,sxu0,sxd,sxd0)
  elif method == "relDiffAsymErrs" : return relDiffAsymErrs(x,x0,sxu,sxu0,sxd,sxd0)
  elif method == "compatAsym"    : return compatAsym(x,x0,sxu,sxu0,sxd,sxd0)
  elif method == "diffPullAsym"  : return diffPullAsym(x,x0,sxu,sxu0,sxd,sxd0)
  else: sys.exit("python/calculate_pulls.py -- Error, pulls method %s not understood!"%method)

def returnTitle(method):
  if   method in ["unconstPullAsym","unconstPull"] : return '(#theta-#theta_{I})/#sigma'
  elif method in ["relDiffAsymErrs"] : return '(#theta-#theta_{I})/#sigma_{I}'
  elif method in ["compat","compatAsym"] : return '#chi^{2} compat.' 
  elif method in ["diffPullAsym","diffPull"] : return '(#theta-#theta_{I})/#sqrt{#sigma_{I}^{2}-#sigma^{2}}' 
  else: sys.exit("python/calculate_pulls.py -- Error, pulls method %s not understood!"%method)

