import pandas as pd
import ktwo19.plotting
import ktwo19.io
from ktwo19.config import bjd0
from matplotlib.pylab import *
import seaborn as sns
sns.set_style('whitegrid')
sns.set_color_codes()
ephem0 = ktwo19.io.load_ephem().copy()
times0 = ktwo19.io.load_times()
times0 = times0.dropna(subset=['inst'])
times0 = times0[times0.inst.str.contains('K2|Spitzer|FLWO|TRAPPIST|MuSCAT')]
times0.index =times0.i_planet

import ttv.plotting

def predict_spitzer_2017_fall(i_planet, plot_djh=False, plot_spitzer_window=False):
    times_djh = ktwo19.io.load_djh()
    pad_start = 2
    pad_stop = 2
    spitzer_start = "2017-08-30 00:00:00"
    spitzer_stop = "2017-10-07 18:39:00"
    if i_planet==1:
        T14 = 3.237
    elif i_planet==2:
        T14 = 3.8

    times = times0.ix[i_planet]
    per = ephem0.ix[i_planet,'per']
    T = ephem0.ix[i_planet,'T']

    i_epoch, tc = ktwo19.io.model_from_samples('ttvfast-npar9',i_planet)
    quants = pd.DataFrame(tc).quantile([0.05,0.15,0.5,0.85,0.95])
    quants = quants.T
    tc_lo = quants[0.05]
    tc_hi = quants[0.95]

    window = ttv.predict.windows(
        i_epoch,tc_lo,tc_hi,T14,spitzer_start,spitzer_stop,pad_start,pad_stop,
        bjd0=bjd0,full_output=True
    )
    window['tmid'] = 0.5*(window['tstart'] + window['tstop'] )
    window['tdur'] = window['tstop'] - window['tstart'] 

    ttv.plotting.plot_omc(i_epoch,tc.T,per,T,color='r',alpha=0.1)
    _ = ttv.plotting.plot_omc(i_epoch,tc_lo,per,T,color='k',label='tc 5%-ile')
    _ = ttv.plotting.plot_omc(i_epoch,tc_hi,per,T,color='b',label='tc 95%-ile')
    ttv.plotting.errorbar_omc(
        times.i_epoch,times.tc,per,T,color='k',fmt='o',yerr=times.tc_err*5,
        label='Observed times\n(5$\sigma$ errorbars)'
    )

    if plot_djh:
        _times = times_djh.ix[i_planet]
        ttv.plotting.errorbar_omc(
            _times.i_epoch,_times.tc,per,T,yerr=_times.tc_err,color='m',
            fmt='.',label='DJH Pred.'
        )

    if plot_spitzer_window:
        _ = ttv.plotting.errorbar_omc(
            window.i_epoch,window.tmid,per,T,yerr=window.tdur/2,color='k',
            fmt='.',ms=0,label='Spitzer Window'
        )

    legend(loc='best')
    print window['tstart_iso tstop_iso tobs_hr tc_lo tc_hi'.split()]
    xlim(2000,3500)
    xlabel('BJD - %i' % bjd0)
    ylabel('TTV (days)')
