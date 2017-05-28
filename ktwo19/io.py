import os
from cStringIO import StringIO

import numpy as np
import pandas as pd

from .config import DATADIR, bjd0, M_earth
import cPickle as pickle
import pymultinest
import re
import ttv.lithwick
from cStringIO import StringIO as sio


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
        outputfiles_basename=u'analysis/mnest-chains/ttvfast-npar9_lp=300/1-'
        an = pymultinest.analyse.Analyzer(
            8, outputfiles_basename=outputfiles_basename
        )
        cols= ['per1', 'ecosw1', 'esinw1', 'tc1', 'mr2', 'per2', 'ecosw2',
       'esinw2', 'tc2']
        samples = pd.DataFrame(
            an.get_equal_weighted_posterior()[:,:-1],columns=cols
        )
        samples = samples.rename(columns={'T1':'tc1','T2':'tc2'})

    if method=='ttvfast-npar9-secosw':
        outputfiles_basename=u'analysis/mnest-chains/ttvfast-npar9-secosw/1-'
        an = pymultinest.analyse.Analyzer(
            8, outputfiles_basename=outputfiles_basename
        )
        cols= ['per1', 'secosw1', 'sesinw1', 'tc1', 'mr2', 'per2', 'secosw2',
       'sesinw2', 'tc2']
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

def load_djh():
    s = """\
    tc               tc_err
    2457898.54880    0.00459 
    2457906.46919    0.00569 
    2457914.39057    0.00674 
    2457922.31341    0.00776 
    2457930.23357    0.00881 
    2457938.15467    0.00981 
    2457946.07734    0.01071 
    2457953.99728    0.01167 
    2457961.91812    0.01261 
    """
    s = pd.read_table(sio(s),sep='\s*')
    s['tc'] -=bjd0
    s['i_epoch'] = s.index
    s['i_epoch'] += 137
    s['i_planet'] = 1
    times_djh = s.copy()


    s = """\
    tc               tc_err
    2457900.04646    0.02629 
    2457911.94202    0.03194 
    2457923.83491    0.03837 
    2457935.73184    0.04370 
    2457947.62592    0.04959 
    2457959.52420    0.05443 
    2457971.41940    0.05959 
    """
    s = pd.read_table(sio(s),sep='\s*')
    s['tc'] -=bjd0
    s['i_epoch'] = s.index
    s['i_epoch'] += 91
    s['i_planet'] = 2

    times_djh = times_djh.append(s)
    times_djh.index = times_djh.i_planet
    return times_djh
