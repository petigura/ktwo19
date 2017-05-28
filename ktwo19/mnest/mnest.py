from __future__ import absolute_import, unicode_literals, print_function
from matplotlib import pylab as plt

import pymultinest
import math
import os
import threading, subprocess
from sys import platform
import ktwo24.io
import ttv.lithwick
import lmfit
import numpy as np
import json

#import ttv.lithwick.nresiduals

def show(filepath): 
   """ open the output (pdf) file for the user """
   if os.name == 'mac' or platform == 'darwin': 
      subprocess.call(('open', filepath))
   elif os.name == 'nt' or platform == 'win32': 
      os.startfile(filepath)
   elif platform.startswith('linux'): 
      subprocess.call(('xdg-open', filepath))

def run_mcmc(myloglike, myprior, nparams, parameters, **kwargs):
   import matplotlib.pyplot as plt
   outputfiles_basename = kwargs['outputfiles_basename']

   nrows = 2
   ncols = int(np.ceil(1.0 * nparams / nrows))
   fig, axL = plt.subplots(nrows=nrows,ncols=ncols,figsize=(20,10),sharey=True)
   axL = axL.flatten()
   for i in range(nparams):
      plt.sca(axL[i])
      plt.title(parameters[i])
      fig.set_tight_layout(True)
      plt.yscale('symlog',ythresh=0.05)
      plt.ylim(-1e10,-1)

   # we want to see some output while it is running
   progress = pymultinest.ProgressPlotter(
      n_params=nparams, outputfiles_basename=outputfiles_basename, axL=axL
   )


   progress.start()
   pdffn = "{}phys_live.points.png".format(outputfiles_basename)
   # delayed opening
   threading.Timer(10, show, [pdffn]).start() 


   # run MultiNest
   pymultinest.run(myloglike, myprior, nparams, **kwargs)


   progress.stop()


   # lets analyse the results
   a = pymultinest.Analyzer(
      n_params=nparams, outputfiles_basename=outputfiles_basename
   )
   s = a.get_stats()

   # store name of parameters, always useful
   with open('%sparams.json' % a.outputfiles_basename, 'w') as f:
      json.dump(parameters, f, indent=2)

   # store derived stats
   with open('%sstats.json' % a.outputfiles_basename, mode='w') as f:
      json.dump(s, f, indent=2)

   print()
   print("-" * 30, 'ANALYSIS', "-" * 30)
   print("Global Evidence:\n\t%.15e +- %.15e" % ( s['nested sampling global log-evidence'], s['nested sampling global log-evidence error'] ))

   plt.clf()
   # Here we will plot all the marginals and whatnot, just to show off
   # You may configure the format of the output here, or in matplotlibrc
   # All pymultinest does is filling in the data of the plot.
   # Copy and edit this file, and play with it.
   p = pymultinest.PlotMarginalModes(a)
   plt.figure(figsize=(5*nparams, 5*nparams))
   #plt.subplots_adjust(wspace=0, hspace=0)
   for i in range(nparams):
      plt.subplot(nparams, nparams, nparams * i + i + 1)
      p.plot_marginal(i, with_ellipses = True, with_points = False, grid_points=50)
      plt.ylabel("Probability")
      plt.xlabel(parameters[i])

      for j in range(i):
         plt.subplot(nparams, nparams, nparams * j + i + 1)
         p.plot_conditional(i, j, with_ellipses = False, with_points = True, grid_points=30)
         plt.xlabel(parameters[i])
         plt.ylabel(parameters[j])

   plt.savefig("chains/marginals_multinest.pdf") #, bbox_inches='tight')
   show("chains/marginals_multinest.pdf")

   for i in range(nparams):
      outfile = '%s-mode-marginal-%d.pdf' % (a.outputfiles_basename,i)
      p.plot_modes_marginal(i, with_ellipses = True, with_points = False)
      plt.ylabel("Probability")
      plt.xlabel(parameters[i])
      plt.savefig(outfile, format='pdf', bbox_inches='tight')
      plt.close()

      outfile = '%s-mode-marginal-cumulative-%d.pdf' % (a.outputfiles_basename,i)
      p.plot_modes_marginal(i, cumulative = True, with_ellipses = True, with_points = False)
      plt.ylabel("Cumulative probability")
      plt.xlabel(parameters[i])
      plt.savefig(outfile, format='pdf', bbox_inches='tight')
      plt.close()

   print("Take a look at the pdf files in chains/") 



 
