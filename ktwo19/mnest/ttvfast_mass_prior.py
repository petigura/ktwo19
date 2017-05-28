import numpy as np
import lmfit
from scipy.optimize import minimize
import pandas as pd
from astropy import constants as c 

import ktwo19.io
import ttv.ttvfastmod
from ktwo19.mnest.priors import UniformPrior, GaussianPrior, PriorCollection
M_earth = (c.M_earth / c.M_sun).value 

class Result(object):
   """Generic picklable bucket"""
   def __init__(self):
      pass

times = ktwo19.io.load_times()
times.index = [times.i_planet,times.i_epoch]
times = times.dropna(subset=['inst'])
times = times[times.inst.str.contains('K2|FLWO|TRAPPIST|MuSCAT|Spitzer')]

parameters_full = [
   "mass1", "per1", "ecosw1", "esinw1", "inc1", "node1", "T1",
   "mass2", "per2", "ecosw2", "esinw2", "inc2", "node2", "T2",
   "stellar_mass"
]
parameters_full = np.array(parameters_full)

fixed = [
   False, False, False, False, True, True, False,
   False, False, False, False, True, True, False,
   True,
]
fixed = np.array(fixed)

fixed_val = [
   np.nan, np.nan, np.nan, np.nan, 90.0, 180.0, np.nan,
   np.nan, np.nan, np.nan, np.nan, 90.0, 180.0, np.nan,
   0.93,
]
fixed_val = np.array(fixed_val)

print "variable pararmeters:"
print parameters_full[~fixed]

print "fixed pararmeters:"
print parameters_full[fixed]

# number of dimensions our problem has
parameters = parameters_full[~fixed]
nparams = len(parameters)
print nparams

dt = 0.3
tstart = 1976  # Chosen so the first transit is the first K2 transit.
tstop = 2400
ttvmod = ttv.ttvfastmod.TTVFastModel(
   2, tstart, dt, tstop, fixed, fixed_val, 
   basis="mass per ecosw esinw inc node tc"
)
ntransits = len(times)
times = np.array(times.sort_values(by='tc').to_records(index=False))

def nresid(params):
   times_mod = ttvmod.compute(params)
   times_mod = ttv.ttvfastmod.parse_ttv_results(times_mod, ttvmod.nplanet)
   times_mod = pd.DataFrame(times_mod)
   times_mod.index = [times_mod.i_planet,times_mod.i_epoch]
   _resid = (times.tc - times_mod.time).dropna()
   _nresid = _resid / times.tc_err.ix[_resid.index]
   _nresid = np.array(_nresid)

   return _nresid

def chisq(params):
   _nresid = nresid(params)
   if len(_nresid) < ntransits:
      return 1e90
   _chisq = np.sum(_nresid**2)
   return _chisq


#      GaussianPrior(20.8858066, 0.001), # per1 
per1 = 7.920994 # narita
per2 = 12.0028 # narita

#per1 = 7.91940 # sinukoff
#per2 = 11.90715 # sinukoff
tc1 = 1980.38403
tc2 = 1984.27227
priors = PriorCollection(   
   [
      UniformPrior(per1 - 0.01, per1 + 0.01), # per1 
      UniformPrior(-0.1,0.1), # ecosw1
      UniformPrior(-0.1,0.1), # esinw1
      GaussianPrior(tc1, 0.01), # tc1 
      UniformPrior(0.1,0.5), # mr2
      UniformPrior(per2 -0.001, per2 + 0.001), # per2
      UniformPrior(-0.1,0.1), # ecosw2
      UniformPrior(-0.1,0.1), # esinw2
      GaussianPrior(tc2, 0.01), # tc2
   ]
)


def priorcube(cube, ndim, nparams):
   """maps the unit cube to the physical parameters how we set bounds"""
   priors.prior_transform(cube, ndim, nparams)

def loglikecube(cube, ndim, nparams):
   # update when changing number of variable parameters
   params = np.array(
      [
         cube[0],
         cube[1],
         cube[2],
         cube[3],
         cube[4],
         cube[5],
         cube[6],
         cube[7],
         cube[8],
      ]
   )
   _loglike = -0.5 * chisq(params)
   return _loglike

def cubemap(u, x0, w):
   # To have the cube stay between x0 - w and x0 + w, perform the
   # following mapping:
   # y = (2*w) * u + x0 - w
   return 2*w*(u - 0.5) + x0

def loglike(params):
   _loglike = -0.5 * chisq(params)
   _loglike+=prior(params)
   return _loglike

res = Result()
res.ttvmod = ttvmod
