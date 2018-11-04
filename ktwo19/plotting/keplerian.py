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

def fig_rv():
    plt.rcParams['lines.markersize'] = 3
    post = ktwo19.keplerian.maximum_a_posteriori()
    plotter = orbit_plots.GPMultipanelPlot(
        post,subtract_orbit_model=False,nobin=True,
        phase_nrows=None, phase_ncols=3,
    )
    fig = plt.figure(figsize=(8,5.5))
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=3)
    ax2 = plt.subplot2grid((2, 3), (1, 0),)
    ax3 = plt.subplot2grid((2, 3), (1, 1),)
    ax4 = plt.subplot2grid((2, 3), (1, 2),)
    axL = [ax1,ax2,ax3,ax4]
    plt.rc('font',size=9)
    plotter.nobin = True
    plotter.legend = False
    plotter.epoch = 2454833
    plotter.phasetext_size = 'small'
    plt.sca(ax1)
    plotter.plot_timeseries()
    radvel.plot.labelfig(97)

    plt.sca(ax2)
    plotter.plot_phasefold(98, 1)
    plt.sca(ax3)

    plotter.plot_phasefold(99, 2)
    plt.sca(ax4)
    plotter.plot_phasefold(100, 3)


    plt.setp(axL[2:],ylabel='')
    _ = plt.setp(axL[1:],ylim=(-23,23),yticks=[-20,-15,-10,-5,0,5,10,15,20])
    _ = plt.setp(axL[0],ylim=(-38,38),xlim=(2200,3300))
    plt.tight_layout(True)
    plt.subplots_adjust(wspace=0.2)

