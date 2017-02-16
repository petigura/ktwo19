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
    times_k2 = []
    csvfn = 'times_pipe=paper_starname=K2-24_i-planet=1.csv'
    csvfn = os.path.join(DATADIR, 'k2/',csvfn)
    times = pd.read_csv(csvfn,index_col=None)
    times['planet'] = 'b'
    times_k2.append(times)
    
    csvfn = 'times_pipe=paper_starname=K2-24_i-planet=2.csv'
    csvfn = os.path.join(DATADIR, 'k2/',csvfn)
    times = pd.read_csv(csvfn,index_col=None)
    times['planet'] = 'c'
    times_k2.append(times)

    times_k2 = pd.concat(times_k2)
    times_k2['inst'] = 'K2'
    times_k2['t0']+=bjd0
    times_k2 = times_k2.rename(
        columns={'itrans':'i_epoch','t0':'tc','t0_err':'tc_err'}
    )

    # Next load spitzer times
    csvfn = 'transit-times-spitzer.csv'
    csvfn = os.path.join(DATADIR, 'spitzer/', csvfn)
    times_spitzer = pd.read_csv(csvfn)
    times_spitzer = pd.read_csv(csvfn)

    times_spitzer.index = times_spitzer.planet
    times_spitzer['inst'] = 'Spitzer'
    times = pd.concat([times_k2,times_spitzer])
    times['i_planet'] = times.planet.replace('b',1).replace('c',2).astype(int)
    times.index = [times.planet,times.i_epoch]

    times.ix[('c',15),'tc']+=(5*42.3633)

    times.index = times.planet
    times['tc'] -= bjd0
    g = times.groupby('planet')
    #per = g.agg(per_linear)['i_epoch']
    #t0 = g.agg(t0_linear)['i_epoch']

    #times['per_linear'] = per
    #times['t0_linear'] = t0
    times.index = [times.planet,times.i_epoch]
    times = times['i_planet i_epoch tc tc_err planet inst notes'.split()]
    return times

def load_ephem(method='lithwick'):
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
        print ephem
    else:
        raise ValueError, "method {} not supported".format(method)

    return ephem

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


def read_dai():
    basedir = os.path.join('/Users/petigura/Research/subsaturn/subsaturn/data/rv/literature/Dai16/')
    df = []
    for starname in 'K2-19 K2-24 K2-32'.split():
        fn = os.path.join(basedir,'{}.txt'.format(starname) )
        temp = pd.read_table(
            fn, sep='&', prefix='x', header=None, 
            names='t mnvel errvel bis tel'.split()
        )
        temp['starname'] = starname
        df.append(temp)
    df = pd.concat(df)
    df.tel = df.tel.replace(0,'pfs').replace(1,'harps')
    return df

def load_samples(method):
    if method=='ttvfast-npar=11':
        samples = pd.read_hdf(
            'analysis/emcee-ttvfast-npar=11_nburn=1000/chain.hdf','chain'
        )        
        samples = samples[(samples.mass1 / M_earth) < 30]

    if method=='lithwick2':
        outputfiles_basename=u'analysis/mnest-chains/lithwick2-lp=1000/1-'
        #outputfiles_basename=u'analysis/mnest-chains/test/1-'
        an = pymultinest.analyse.Analyzer(
            8, outputfiles_basename=outputfiles_basename
        )
        cols= 'per1 tc1 mu1 per2 tc2 mu2 ReZfree ImZfree'.split()
        samples = pd.DataFrame(
            an.get_equal_weighted_posterior()[:,:-1],columns=cols
        )
        samples = samples.rename(columns={'T1':'tc1','T2':'tc2'})

    if method=='lithwick-bootstrap':
        samples = pd.read_hdf('analysis/lithwick-bootstrap.hdf','samples')
        samples = samples.rename(columns={'T1':'tc1','T2':'tc2'})
        samples = samples['per1 tc1 ReV1 ImV1 per2 tc2 ReV2 ImV2'.split()]
        
    if method=='lithwick2-muprior':
        outputfiles_basename=u'analysis/mnest-chains/lithwick2-muprior-lp=1000/1-'
        an = pymultinest.analyse.Analyzer(
            8, outputfiles_basename=outputfiles_basename
        )
        cols= 'per1 tc1 mu1 per2 tc2 mu2 ReZfree ImZfree'.split()
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

    if method=="lithwick2":
        import ktwo24.mnest.mnest_lithwick2 as mod
        ttvmod = mod.ttvmod
        ttvmod.basis = 'per1 tc1 mu1 per2 tc2 mu2 ReZfree ImZfree'
        i_epoch = np.arange(200)
        i_planet = np.ones_like(i_epoch) * i_planet
        def model(i_samp):
            tmod = ttvmod.compute(samp[i_samp], i_planet, i_epoch)
            tmod = pd.DataFrame(tmod, columns=['tc'])
            tmod['i_epoch'] = i_epoch
            tmod.index = tmod.i_epoch
            return np.array(tmod.tc)

    if method=="ttvfast-npar=11":
        import ktwo24.mnest.emcee_ttvfast_npar11 as mod
        ttvmod = mod.ttvmod
        ttvmod.tstart = 2065
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


    return i_epoch, tc
