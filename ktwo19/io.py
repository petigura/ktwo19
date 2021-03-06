import sys
import os
from cStringIO import StringIO as sio

import numpy as np
import pandas as pd
from scipy import optimize
import scipy.stats
import cpsutils.io
from astropy.time import Time
from astropy import constants as c
from astropy import units as u
from matplotlib.pylab import *
sys.path.append('/Users/petigura/Research/subsaturn/')
import subsaturn.lopez

import ktwo19.phodymm
import ktwo19.keplerian 
#import ktwo19.photometry
#import radvel.utils
from .config import bjd0

DATADIR = os.path.join(os.path.dirname(__file__),'../data/')
def load_table(table, cache=0, cachefn='load_table_cache.hdf', verbose=False):
    """Load tables

    Args:
        table (str): name of table. must be one of


        cache (Optional[int]): whether or not to use the cache
            - 0: don't use the cache recreate all files
            - 1: read from cache
            - 2: write tables to cache

    Returns:
        pandas.DataFrame: table

    """
    if cache==1:
        try:
            df = pd.read_hdf(cachefn,table, mode='a')
            print "read table {} from {}".format(table,cachefn)
            return df

        except IOError:
            print "Could not find cache file: %s" % cachefn
            print "Building cache..."
            cache=2
    
        except KeyError:
            print "Cache not built for table: %s" % table
            print "Building cache..."
            cache=2

    if cache==2:
        df = load_table(table, cache=False)
        print "writing table {} to cache".format(table)
        df.to_hdf(cachefn,table)
        return df

    elif table=='stellar':
        fn = os.path.join(DATADIR, 'data.xlsx')
        df = pd.read_excel(fn,'stellar',squeeze=True,header=None,index_col=0,usecols=1)

    elif table=='ephem-sinukoff16':
        fn = 'data.xlsx'
        fn = os.path.join(DATADIR, fn)
        df = pd.read_excel(fn,sheet_name=table,usecols=[0,1],index_col=0,squeeze=True)
        for i in range(1,4):
            df['tc%i' % i] += (2456000 - 2454833) # Kepler epoch

    elif table=='ktwo-everest':
        df = ktwo19.photometry._everest()

    elif table=='rv-full':
        vstfn = os.path.join(DATADIR,'rv/vstepic201505350.dat')
        vst, meta = cpsutils.io.read_vst(vstfn,full_output=True, verbose=0)
        sval = cpsutils.io.loadsval()
        sval['obs'] = sval.obs.str.replace('bj','rj')
        sval = sval.rename(columns={'obs':'obnm'})
        vst = pd.merge(vst,sval['obnm sval'.split()])

        vst = vst['jd day mnvel errvel sval obnm cts'.split()]
        vst = vst.rename(columns={'jd':'time'})
        vst['tel'] = 'j'
        vst['time'] -= bjd0
        df = vst

    elif table=='rv':
        df = load_table('rv-full')
        omit=['rj197.149'] # 2 deg from 90% full moon, check sky in obs
        # omit=['rj197.345'] # 14 deg from 80% full moon 
        # omit=['rj196.319'] # 22 deg from full moon, possible light cirrus, 50k to be safe 
        # omit = ['rj196.319','rj197.149','rj197.345', 'rj200.88']#,'rj218.227','rj219.300','rj221.153','rj222.136','rj222.327','rj224.149']
        for i in omit:
            if i in df.obnm.tolist():
                df=df.drop(df.index[df.obnm==i]).reset_index(drop=True) 
                print 'excluding obs: ' + i


    elif table=='rv-trend-removed':
        P, post = radvel.utils.initialize_posterior('analysis/radvel/ktwo19_npl=2-cc.py')
        res  = optimize.minimize(post.neglogprob_array, post.get_vary_params(), method='Nelder-Mead',) # MAP fit

        print "MAP fit"
        print post
        post.params['dvdt'].value = 0 
        post.params['gamma_j'].value = 0 

        print "Zeroing out dvdt and trend"
        print post
        model = post.likelihood.model
        rvtimes = post.likelihood.x # times where model is computed  
        rvmod = model(rvtimes) # 
        rawresid = post.likelihood.residuals() # array
        rv_dt = rawresid + rvmod # residuals plus model without trend
        rv_dt_err = post.likelihood.yerr
        
        df = dict(time=rvtimes,mnvel=rv_dt,errvel=rv_dt_err)
        df = pd.DataFrame(df)
        df2 = load_table('rv')
        df = pd.merge(df,df2['obnm tel sval time'.split()])
        df = df['tel obnm time mnvel errvel sval'.split()]

    elif table=='stellar':
        fn = os.path.join(DATADIR, 'data.xlsx')
        df = pd.read_excel(
            fn,'stellar-sme',squeeze=True,header=None,index_col=0,usecols=[0,1]
        )

    elif table=='times':
        fn = os.path.join(DATADIR, 'data.xlsx')
        df = pd.read_excel(fn,'transit-times')
        df = df[df.include==1]

        jd = df.tc
        df['day'] = Time(df.tc,format='jd',out_subfmt='date').iso
        df['tc'] -= 2454833


    elif table=='phot-gp':
        fn = os.path.join(DATADIR, 'data.xlsx')
        df = pd.read_excel(
            fn,sheetname='phot-gp',squeeze=True,header=None,index_col=0
        )

    elif table=='keplerian-samples':
        df = ktwo19.keplerian.mcmc()

    elif table=='keplerian-samples-derived':
        star = load_table('stellar')
        samp = load_table('keplerian-samples',cache=1)
        ephem = load_table('ephem-sinukoff16')

        size = len(samp)
        smass = np.random.normal(loc=star.smass,scale=star.smass_err,size=size)
        srad = np.random.normal(loc=star.srad,scale=star.srad_err,size=size)
        samp['smass'] = smass
        for i in range(1,4):
            si = str(i)
            K = samp['k%i' % i]
            per = ephem['per%i'%i]
            masse = radvel.utils.Msini(K, per, smass, 0, Msini_units='earth')
            mu = np.array((masse * c.M_earth) / (smass *c.M_sun))
            samp['mpsini%i' %i] = masse
            samp['musini%i' %i] = mu * 1e6

            fn = '201505350.0{}-mcmc-samples.csv.gz'.format(i)
            fn = os.path.join('data/livingston-lc-fits/',fn)
            lcsamp = pd.read_csv(fn)            
            ror = lcsamp.sample(size,replace=True)['k'] 
            samp['prad%i' % i] = np.array(ror * (srad * c.R_sun) / c.R_earth)

            '''
            pmass = df['masse%i' % i]
            df['prad'+si] = prad
            rho = (np.array(pmass)*c.M_earth) / (4.0/3.0 * np.pi * (np.array(prad) * c.R_earth)**3)
            rho = np.array(rho.to(u.g * u.cm**-3))
            df['rho'+si] = rho

            '''
        df = samp

    elif table=='photodyn-samples':
        basedir = 'analysis/photodyn/runs/'

        # Use for purely photdynamical fit
        '''
        demcmcfname = 'K2-19_e-uniform_Omega-vary_no-RV/demcmc_k2-19_massprior.out'
        fname = 'K2-19_e-uniform_Omega-vary_no-RV/k2-19.in'
        '''


        demcmcfname  = 'K2-19_e-uniform_Omega-vary/demcmc_k2-19_massprior.out'
        fname = 'K2-19_e-uniform_Omega-vary/k2-19.in'

        demcmcfname = os.path.join(basedir, demcmcfname)
        fname = os.path.join(basedir, fname)
        nburn = 10000 
        df = ktwo19.phodymm.read_phodymm_format(fname, demcmcfname, nburn)

    elif table=='fenv-samples':
        lopi = subsaturn.lopez.LopezInterpolator()
        df = load_table('photodyn-samples')

        stellar = load_table('stellar')
        nsamp = 1000
        teff = scipy.stats.norm.rvs(stellar.steff, stellar.steff_err,nsamp)
        smass = scipy.stats.norm.rvs(stellar.smass, stellar.smass_err,nsamp)
        srad = scipy.stats.norm.rvs(stellar.srad, stellar.srad_err,nsamp)
        test = pd.DataFrame(index=range(nsamp))

        for i in [1,2,3]:
            pmass = df['masse%i' % i]
            prad = df['prad%i' %i ]
            per = np.array(df['per%i' %i].sample(nsamp))
            a = per2a(per*u.day, smass*u.M_sun ).to(u.au).value
            _teq = teq(teff, srad, a)

            plnt = dict(
                pl_masse=np.median(pmass),
                pl_masseerr1=np.std(pmass),
                pl_masseerr2=-np.std(pmass),
                pl_rade=np.median(prad),
                pl_radeerr1=np.std(prad),
                pl_radeerr2=-np.std(prad),
                pl_teq=np.median(_teq),
                pl_teqerr1=np.std(_teq),
                pl_teqerr2=-np.std(_teq),
                age=5
            ) 
            plnt = pd.Series(plnt)
            temp = subsaturn.lopez.sample_ss_cmf(lopi,plnt,age=5,size=nsamp)
            test['fenv%i' % i] = 1 - temp['cmf']
            test['teq%i' % i ] = _teq
            for k in 'mcore menv cmf'.split():
                test[k+str(i)] = temp[k]

            for k in 'fenv cmf'.split():
                test[k+str(i)] *= 100
        df = test



    else:
        assert False, "Table {} not valid table name".format(table)

    return df


def teq(teff, rstar, a):
    """
    teff : effective temperature (K)
    rstar: stellar radii (solar radii)
    a: semi-major axis AU
    """
    teq_earth = 278
    _teq = teq_earth * (teff / 5770.0) * (rstar)**0.5 * a**-0.5
    return _teq

def per2a(per, mstar):
    a = (per**2 / 4 / np.pi**2 * c.G * mstar)**(1/3.)
    return a

def planet_line(line,pnum):
    line = line.split('\t')
    d = {}
    d['per%i' % pnum] = line[1] # period (days)
    d['tc%i' % pnum] = line[2] # time of conjunction (days)
    d['secosw%i' % pnum] = line[3] # sqrt(e)cos(w) w = argument of peri
    d['sesinw%i' % pnum] = line[4] # sqrt(e)cos(w)
    d['inc%i' % pnum] = line[5] # inclination deg
    d['Omega%i' % pnum] = line[6] # Longitude of ascending node (deg)
    d['massj%i' % pnum] = line[7] # Mass in jupiter units
    d['rrat%i' % pnum] = line[8][:-2] # Radius ratio
    return d

def load_ephem(method='linear'):
    fn = os.path.join(DATADIR, 'data.xlsx')
    df = pd.read_excel(fn,'transit-times')
    df = df[(df.include==1) | (df.inst=='K2')]
    df['tc'] -= 2454833
    times = df
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
