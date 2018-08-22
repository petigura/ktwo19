import numpy as np
import pandas as pd
import os
import radvel
from radvel.plot import orbit_plots, mcmc_plots
from scipy import optimize
import ktwo19
import ktwo19.io
import ktwo19.keplerian
from matplotlib import pylab as plt

def rv():
    post = ktwo19.keplerian.maximum_a_posteriori()
    plotter = orbit_plots.GPMultipanelPlot(
        post,subtract_orbit_model=False,nobin=True
    )
    plotter.plot_multipanel()
    fig = plt.gcf()
    axL = fig.get_axes()
    plt.setp(axL[0:2],ylim=(-30,30))
