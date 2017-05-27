"""
Example:


"""

import os
from sys import platform
import ktwo19.io
import ktwo19.mnest.mnest
from matplotlib import pylab as plt
from argparse import ArgumentParser
import cPickle as pickle

fig = plt.figure(figsize=(10,10))
fig.set_tight_layout(True)

def main():
    psr = ArgumentParser(description='')
    psr.add_argument('method')
    psr.add_argument('id')
    psr.add_argument('--livepoints', default=1000, type=int)
    psr.add_argument('--ins', default=False, action='store_true')
    psr.add_argument('--ceff', default=False, action='store_true')
    psr.add_argument('--etol', default=5,type=float)
    args = psr.parse_args()

    outputfiles_basename = 'analysis/mnest-chains/{}/1-'.format(args.id)
    dirname = os.path.dirname(outputfiles_basename)
    print "mnest output files: {}".format(dirname) 
    if not os.path.exists(dirname): 
       os.makedirs(dirname)
        
    if args.method=="ttvfast-mass-prior":
        import ktwo19.mnest.mnest_ttvfast_mass_prior as mod
    if args.method=="ttvfast-npar9":
        import ktwo19.mnest.ttvfast_npar9 as mod
    else:
        assert False, "method not supported"


    ktwo19.mnest.mnest.run_mcmc(
        mod.loglikecube, mod.priorcube, mod.nparams, mod.parameters, 
        importance_nested_sampling = args.ins, 
        resume=True, verbose=True, sampling_efficiency='model', 
        n_live_points=args.livepoints, 
        outputfiles_basename=outputfiles_basename
    )

if __name__=="__main__":
    main()


