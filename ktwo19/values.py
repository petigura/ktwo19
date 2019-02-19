import imp
from collections import OrderedDict
import glob
import os 

import pandas as pd
import numpy as np
import radvel.driver
import radvel.utils
from astropy import constants as c

import ktwo19.io
import cpsutils.io

def val_stat(return_dict=False):
    d = OrderedDict()
    vst = ktwo19.io.load_table('rv')
    d['rv-n'] = len(vst)
    d['rv-start'] = vst.iloc[0]['day']
    d['rv-stop'] = vst.iloc[-1]['day']
    d['rv-errvel-min'] = "{:.1f}".format(vst.errvel.min())
    d['rv-errvel-max'] = "{:.1f}".format(vst.errvel.max())
    d['rv-errvel-med'] = "{:.1f}".format(vst.errvel.median())
    d['rv-snr-min'] = "{:.0f}".format(np.sqrt(vst.cts.min()))
    d['rv-snr-max'] = "{:.0f}".format(np.sqrt(vst.cts.max()))
    d['rv-snr-med'] = "{:.0f}".format(np.sqrt(vst.cts.median()))

    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d

    return lines

def val_sys(return_dict=False):
    d = OrderedDict()
    df = ktwo19.io.load_table('stellar',cache=1)
    for k in df.index:
        d[k] = df[k]
    d['steff'] = "{:.0f}".format(df.ix['steff'])
    d['steff_err'] = "{:.0f}".format(df.ix['steff_err'])

    df = ktwo19.io.load_table('phot-gp',cache=1)
    d['phot-gp-eta1'] = "{:.1f}".format(df.ix['eta1'] * 1e2)
    d['phot-gp-eta1_err'] = "{:.1f}".format(df.ix['eta1_err'] * 1e2)
    d['phot-gp-eta2'] = "{:.0f}".format(df.ix['eta2'])
    d['phot-gp-eta2_err'] = "{:.0f}".format(df.ix['eta2_err'])
    d['phot-gp-eta3'] = "{:.1f}".format(df.ix['eta3'])
    d['phot-gp-eta3_err'] = "{:.1f}".format(df.ix['eta3_err'])
    d['phot-gp-eta4'] = "{:.2f}".format(df.ix['eta4'])
    d['phot-gp-eta4_err'] = "{:.2f}".format(df.ix['eta4_err'])

    fmt = OrderedDict()
    post = ktwo19.keplerian.posterior()
    d['rv-per1'] = "{:.5f}".format(post.params['per1'].value)
    d['rv-per2'] = "{:.5f}".format(post.params['per2'].value)
    d['rv-per3'] = "{:.5f}".format(post.params['per3'].value)
    d['rv-tc1'] = "{:.3f}".format(post.params['tc1'].value)
    d['rv-tc2'] = "{:.3f}".format(post.params['tc2'].value)
    d['rv-tc3'] = "{:.3f}".format(post.params['tc3'].value)

    chain = ktwo19.io.load_table('keplerian-samples-derived',cache=1,cachefn='load_table_cache-rv.hdf')

    fmt['k1'] = "{:.1f}"
    fmt['k2'] = "{:.1f}"
    fmt['k3'] = "{:.1f}"
    fmt['mpsini1'] = "{:.0f}"
    fmt['mpsini2'] = "{:.0f}"
    fmt['mpsini3'] = "{:.0f}"
    fmt['musini1'] = "{:.0f}"
    fmt['musini2'] = "{:.0f}"
    fmt['musini3'] = "{:.1f}"


    '''
    fmt['prad1'] = "{:.1f}"
    fmt['prad2'] = "{:.1f}"
    fmt['prad3'] = "{:.1f}"
    '''

    fmt['gamma_j'] = '{:.0f}'
    fmt['dvdt'] = '{:.1f}'
    fmt['gp_amp'] = "{:.1f}"
    fmt['gp_explength'] = "{:.1f}"
    fmt['gp_per'] = "{:.1f}"
    fmt['gp_perlength'] = "{:.2f}"

    pre = 'rv-'
    insert_chain_dict(chain, d, fmt, pre=pre) 
    d[pre+'mpsini2_p95']  = "{:.1f}".format(chain['mpsini2'].quantile(0.95))
    d[pre+'mpsini3_p95']  = "{:.1f}".format(chain['mpsini3'].quantile(0.95))

    chain = ktwo19.io.load_table('photodyn-samples',cache=2)
    chain['delta'] = chain.eval('per3/per2 * (2/3) - 1')
    chain['masse3onmasse2'] = chain.eval('masse3/masse2')
    fmt = OrderedDict()
    fmt['mstar'] = "{:.2f}"
    fmt['rstar'] = "{:.2f}"
    fmt['masse1'] = "{:.1f}"
    fmt['masse2'] = "{:.0f}"
    fmt['masse3'] = "{:.1f}"
    fmt['per1'] = "{:.4f}"
    fmt['per2'] = "{:.4f}"
    fmt['per3'] = "{:.4f}"
    fmt['tc1'] = "{:.4f}"
    fmt['tc2'] = "{:.4f}"
    fmt['tc3'] = "{:.4f}"
    fmt['secosw1'] = "{:.2f}"
    fmt['secosw2'] = "{:.2f}"
    fmt['secosw3'] = "{:.2f}"
    fmt['sesinw1'] = "{:.2f}"
    fmt['sesinw2'] = "{:.2f}"
    fmt['sesinw3'] = "{:.2f}"
    fmt['inc1'] = "{:.1f}"
    fmt['inc2'] = "{:.1f}"
    fmt['inc3'] = "{:.1f}"
    fmt['Omega1'] = "{:.1f}"
    fmt['Omega2'] = "{:.1f}"
    fmt['Omega3'] = "{:.1f}"
    fmt['omegadeg1'] = "{:.0f}"
    fmt['omegadeg2'] = "{:.0f}"
    fmt['omegadeg3'] = "{:.0f}"
    fmt['omegadiffdeg'] = "{:.0f}"
    fmt['e1'] = "{:.2f}"
    fmt['e2'] = "{:.2f}"
    fmt['e3'] = "{:.2f}"
    fmt['ror1'] = "{:.4f}"
    fmt['ror2'] = "{:.4f}"
    fmt['ror3'] = "{:.4f}"
    fmt['c1'] = "{:.1f}"
    fmt['c2'] = "{:.1f}"
    fmt['prad1'] = "{:.2f}"
    fmt['prad2'] = "{:.1f}"
    fmt['prad3'] = "{:.1f}"

    insert_chain_dict(chain, d, fmt, pre='pd-') 

    chain = ktwo19.io.load_table('fenv-samples',cache=1)
    fmt = OrderedDict()
    fmt['fenv2'] = "{:.0f}"
    fmt['mcore2'] = "{:.0f}"
    fmt['fenv3'] = "{:.0f}"
    fmt['mcore3'] = "{:.1f}"
    fmt['teq1'] = "{:.0f}"
    fmt['teq2'] = "{:.0f}"
    fmt['teq3'] = "{:.0f}"
    
    insert_chain_dict(chain, d, fmt, pre='lopez-') 


    '''
    pre = 'ccc-'
    chain['dvdt'] *= 365
    chain['musini1'] *= 1e6
    chain['musini2'] *= 1e6
    chain['musini3'] *= 1e6

    star = ktwo19.io.load_table('stellar')
    smass = np.random.normal(
        loc=star.smass, scale=star.smass_err, size=len(chain)
    )
    
    chain['a1'] = radvel.utils.semi_major_axis(chain.per1, smass)
    chain['a2'] = radvel.utils.semi_major_axis(chain.per2, smass)
    chain['a3'] = radvel.utils.semi_major_axis(chain.per3, smass)
    insert_chain_dict(chain, d, fmt, pre='ccc-') 

    ## Eccentric model
    chain =  ktwo19.io.load_table('ktwo19_npl=3-eec',cache=1)
    pre = 'eec-'
    chain['e1'] = chain.eval('secosw1**2 + sesinw1**2')
    chain['e2'] = chain.eval('secosw2**2 + sesinw2**2')
    d[pre+'e1_p90']  = "{:.2f}".format(chain['e1'].quantile(0.9))
    d[pre+'e2_p90']  = "{:.2f}".format(chain['e2'].quantile(0.9))
    '''

    if return_dict:
        return d
    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    return lines

def val_keplerian():
    """
    Print values associated with the Lithwick fits
    """
    chain = ktwo19.io.load_table('keplerian-samples',cache=1)
    chain['muppm1'] = chain.mu1 * 1e6
    chain['muppm2'] = chain.mu2 * 1e6
    chain['zmag']= chain.eval('sqrt( rezfree**2 + imzfree**2)') 
    chain['mr2']= chain.eval('mu2 / mu1') 

    fmt = OrderedDict()
    fmt['per1'] = "{:.5f}"
    fmt['tc1'] = "{:.4f}"
    fmt['muppm1'] = "{:.1f}"
    fmt['per2'] = "{:.4f}"
    fmt['tc2'] = "{:.4f}"
    fmt['muppm2'] = "{:.1f}"
    fmt['rezfree'] = "{:.3f}"
    fmt['imzfree'] = "{:.3f}"
    fmt['zmag'] = "{:.3f}"
    fmt['mr2'] = "{:.2f}"

    d = OrderedDict()
    insert_chain_dict(chain, d, fmt)
    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)
    return lines


def insert_chain_dict(chains, d, fmt, pre=''):
    thresh = 1.5
    # If err1 and err2 are consistent to these values, then report
    # symetric errorbars

    for k in fmt.keys():
        keydict = '{}{}'.format(pre,k)
        keychain = '{}'.format(k)
        s = fmt[k]
        chain = chains[keychain]
        q = chain.quantile([0.16,0.50,0.84])

        val = q.ix[0.50]
        err1 = q.ix[0.84] - q.ix[0.50]
        err2 = q.ix[0.16] - q.ix[0.50]
        is_asymetric = (-err2/err1 > thresh) or (-err1/err2 > thresh)
        
        if not is_asymetric:
            err = 0.5 * (err1 - err2)
            err1 = err 
            err2 = -err

        val = s.format(val)
        err1 = s.format(err1)
        err2 = s.format(err2)

        if is_asymetric:
            s = "$%s^{+%s}_{%s}$" % (val, err1, err2)
        else:
            s = "$%s \pm %s$" % (val, err1)

        s = s.replace('++','+')
        d[keydict] = val
        d[keydict+'_err1'] = err1
        d[keydict+'_err2'] = err2
        d[keydict+'_fmt'] = s

