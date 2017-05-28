"""
TTVFast where initial values are determined by minimization

"""
import numpy as np
import lmfit
from scipy.optimize import minimize
import pandas as pd
from astropy import constants as c 
from ktwo19.mnest.priors import UniformPrior, GaussianPrior, PriorCollection

import ktwo19.io
import ttv.ttvfastmod
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
   "mass1", "per1", "secosw1", "sesinw1", "inc1", "node1", "tc1",
   "mr2", "per2", "secosw2", "sesinw2", "inc2", "node2", "tc2",
   "stellar_mass"
]
parameters_full = np.array(parameters_full)

fixed = [
   True, False, False, False, True, True, False,
   False, False, False, False, True, True, False,
   True,
]
fixed = np.array(fixed)

fixed_val = [
   30.0*M_earth, np.nan, np.nan, np.nan, 90.0, 180.0, np.nan,
   np.nan, np.nan, np.nan, np.nan, 90.0, 180.0, np.nan,
   0.93,
]
fixed_val = np.array(fixed_val)

dt = 0.3
tstart = 1976  # Chosen so the first transit is the first K2 transit.
tstop = 3500
assert tstop > times.tc.max(), "Must integrate to beyond last data point"

ttvmod = ttv.ttvfastmod.TTVFastModel(
   2, tstart, dt, tstop, fixed, fixed_val, 
   basis="mr per secosw sesinw inc node tc"
)
ntransits = len(times)
times = np.array(times.sort_values(by='tc').to_records(index=False))

def nresid(params):
   res = ttvmod.compute(params)
   times_mod, rv = ttv.ttvfastmod.parse_ttv_results(res, ttvmod.nplanet)
   i_planet_data = times['i_planet']
   i_epoch_data = times['i_epoch']
   tc_data = times['tc']
   tc_err_data = times['tc_err']
   i_planet_mod = times_mod['i_planet']
   i_epoch_mod = times_mod['i_epoch']
   tc_mod = times_mod['tc']
   _nresid = ttv.ttvfastmod.ttv_residuals(
      i_planet_data, i_epoch_data, tc_data, tc_err_data,
      i_planet_mod, i_epoch_mod, tc_mod, ttvmod.nplanet
   )
   return _nresid

def chisq(params):
   _nresid = nresid(params)
   if len(_nresid) < ntransits:
      return 1e90

   _chisq = np.sum(_nresid**2)
   return _chisq

# number of dimensions our problem has
parameters = parameters_full[~fixed]
nparams = len(parameters)

#      GaussianPrior(20.8858066, 0.001), # per1 
#per1 = 7.920994 # narita
#per2 = 12.0028 # narita

per1 = 7.91940 # sinukoff
per2 = 11.90715 # sinukoff
tc1 = 1980.38403
tc2 = 1984.27227

# ecosw, esinw
p_restrict0 = [
   per1, 0.001, 0.001, tc1,
   0.5, per2, 0.001, 0.001, tc2
]
p_restrict0 = np.array(p_restrict0)
opt_results = minimize(chisq, p_restrict0, method='Nelder-Mead')
p_restrict1 = opt_results.x
print "guess parameters"
print p_restrict0
print chisq(p_restrict0)
print "final parameters"
print p_restrict1
print chisq(p_restrict1)

#per1 = p_restrict1[0]
#tc1 = p_restrict1[3]
#per2 = p_restrict1[5]
#tc2 = p_restrict1[8]


priors = PriorCollection(   
   [
      UniformPrior(per1 - 0.003, per1 + 0.003), # per1 
      UniformPrior(-0.3,0.3), # secosw1
      UniformPrior(-0.3,0.3), # sesinw1
      GaussianPrior(tc1, 0.003), # tc1 
      UniformPrior(0.1,0.5), # mr2
      UniformPrior(per2 , per2 + 0.01), # per2
      UniformPrior(-0.3,0.3), # secosw2
      UniformPrior(-0.3,0.3), # sesinw2
      GaussianPrior(tc2, 0.003), # tc2
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
