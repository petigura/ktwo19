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

def val_sys():
    d = OrderedDict()
    df = ktwo19.io.load_table('stellar',cache=1)
    for k in df.index:
        d[k] = df[k]
    d['steff'] = "{:.0f}".format(df.ix['steff'])
    d['steff_err'] = "{:.0f}".format(df.ix['steff_err'])

    df = ktwo19.io.load_table('phot-gp',cache=1)
    d['phot-gp-eta1'] = "{:.2f}".format(df.ix['eta1'] * 1e2)
    d['phot-gp-eta1_err'] = "{:.2f}".format(df.ix['eta1_err'] * 1e2)
    d['phot-gp-eta2'] = "{:.0f}".format(df.ix['eta2'])
    d['phot-gp-eta2_err'] = "{:.0f}".format(df.ix['eta2_err'])
    d['phot-gp-eta3'] = "{:.0f}".format(df.ix['eta3'])
    d['phot-gp-eta3_err'] = "{:.0f}".format(df.ix['eta3_err'])
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

    chain = ktwo19.io.load_table('keplerian-samples-derived',cache=2)

    fmt['k1'] = "{:.1f}"
    fmt['k2'] = "{:.1f}"
    fmt['k3'] = "{:.1f}"
    fmt['mpsini1'] = "{:.0f}"
    fmt['mpsini2'] = "{:.0f}"
    fmt['mpsini3'] = "{:.0f}"
    fmt['musini1'] = "{:.0f}"
    fmt['musini2'] = "{:.0f}"
    fmt['musini3'] = "{:.1f}"
    fmt['prad1'] = "{:.1f}"
    fmt['prad2'] = "{:.1f}"
    fmt['prad3'] = "{:.1f}"
    fmt['gamma_j'] = '{:.0f}'
    fmt['dvdt'] = '{:.1f}'
    fmt['gp_amp'] = "{:.1f}"
    fmt['gp_explength'] = "{:.0f}"
    fmt['gp_per'] = "{:.1f}"
    fmt['gp_perlength'] = "{:.2f}"
    insert_chain_dict(chain, d, fmt, pre='rv-') 

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
    for k in fmt.keys():
        keydict = '{}{}'.format(pre,k)
        keychain = '{}'.format(k)
        s = fmt[k]
        chain = chains[keychain]
        q = chain.quantile([0.16,0.50,0.84])
        val = s.format(q.ix[0.50])
        err1 = s.format(q.ix[0.84] - q.ix[0.50])
        err2 = s.format(q.ix[0.16] - q.ix[0.50])

        d[keydict] = val
        d[keydict+'_err1'] = err1
        d[keydict+'_err2'] = err2

        s = "$%s^{+%s}_{%s}$" % (val, err1,err2)
        s = s.replace('++','+')
        d[keydict+'_fmt'] = s

