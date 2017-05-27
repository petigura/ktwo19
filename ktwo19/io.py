import os
from cStringIO import StringIO

import numpy as np
import pandas as pd

from .config import DATADIR, bjd0, M_earth
import cPickle as pickle
import pymultinest
import re
import ttv.lithwick

def load_times():
    # First load K2 times
    times = pd.read_excel('data/transit-times.xlsx',sheet='Sheet1')
    times['notes'] = times.notes.fillna('')
    times['tc'] -= bjd0
    return times

def load_mnest_samples(basedir, nparams, parameters):
    """
    Load results from MultiNest run

    Args:
        basedir (str): path to multinest results
        mod (python module):

    Returns:
        samples (pandas.DataFrame): 

    """
    outputfile_basename = os.path.join(basedir,'1-')

    an = pymultinest.analyse.Analyzer(
        nparams, outputfiles_basename=outputfile_basename
    )
    samples = an.get_equal_weighted_posterior()[:,:-1]
    samples = pd.DataFrame(samples,columns=parameters)
    return samples


def load_ephem(method='linear'):
    times = load_times()
    if method=='lithwick':
        res = ttv.lithwick.fit_lithwick(times,2) 
        ephem = pd.DataFrame(index=[1,2])
        ephem['per'] = [res['per1'],res['per2']]
        ephem['T'] = [res['T1'],res['T2']]
        ephem['per'] = ephem['per'].round(6)
        ephem['T'] = ephem['T'].round(3)
        
    elif method=='linear':
        g = times.groupby('i_planet')
        ephem = pd.DataFrame(index=g.first().index)
        ephem['per'] = g.apply(lambda x : np.polyfit(x['i_epoch'],x['tc'],1)[0] )
        ephem['T'] = g.apply(lambda x : np.polyfit(x['i_epoch'],x['tc'],1)[1] )
    else:
        raise ValueError, "method {} not supported".format(method)

    return ephem



def load_samples(method):
    if method=='ttvfast-npar9':
        outputfiles_basename=u'analysis/mnest-chains/ttvfast-npar9/1-'
        an = pymultinest.analyse.Analyzer(
            8, outputfiles_basename=outputfiles_basename
        )
        cols= ['per1', 'ecosw1', 'esinw1', 'tc1', 'mr2', 'per2', 'ecosw2',
       'esinw2', 'tc2']
        samples = pd.DataFrame(
            an.get_equal_weighted_posterior()[:,:-1],columns=cols
        )
        samples = samples.rename(columns={'T1':'tc1','T2':'tc2'})

    
    print "number of samples {}".format(len(samples))
    return samples

def model_from_samples(method, i_planet, nsamp=100):
    samples = load_samples(method)
    samp = np.array(samples.sample(nsamp,random_state=0))
    tcL = []

    if method=="ttvfast-npar9":
        import ktwo19.mnest.ttvfast_npar9 as mod
        ttvmod = mod.ttvmod
        ttvmod.tstart = 1976
        ttvmod.tstop = ttvmod.tstart + 10 * 365

        def _model(i_samp):
            res = ttvmod.compute(samp[i_samp])
            positions, rv = ttv.ttvfastmod.parse_ttv_results(
                res, ttvmod.nplanet
            )
            times_mod = pd.DataFrame(positions)
            times_mod.index = [times_mod.i_planet,times_mod.i_epoch]    
            tmod = times_mod.ix[i_planet]
            return tmod.i_epoch, np.array(tmod.tc)

        i_epoch = _model(0)[0]
        def model(i_samp):
            return _model(i_samp)[1]

    tc = map(model, range(nsamp))
    minlen = min(map(len, tc))
    tc = [_tc[:minlen] for _tc in tc ]
    tc = np.vstack(tc)
    i_epoch = i_epoch[:minlen]
    return i_epoch, tc
