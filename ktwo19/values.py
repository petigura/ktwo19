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
    # load up data from p16
    p16 = ktwo19.io.load_table('petigura16')
    d['p16-ror1'] = p16.ix['pl_ror',1]
    d['p16-ror1_err'] = p16.ix['pl_ror_err',1]
    d['p16-ror2'] = p16.ix['pl_ror',2]
    d['p16-ror2_err'] = p16.ix['pl_ror_err',2]

    vst, meta = cpsutils.io.read_vst('data/rv/vstepic203771098.dat',full_output=True, verbose=0)
    d['rv-n'] = len(vst)
    d['rv-start'] = vst.iloc[0]['day']
    d['rv-stop'] = vst.iloc[-1]['day']
    d['rv-errvel-min'] = "{:.1f}".format(vst.errvel.min())
    d['rv-errvel-max'] = "{:.1f}".format(vst.errvel.max())

    ephem = ktwo19.io.load_table('ephem-lithwick')
    d['ref-ephem-per-1'] = ephem.ix[1,'per']
    d['ref-ephem-tc-1'] = ephem.ix[1,'T']
    d['ref-ephem-per-2'] = ephem.ix[2,'per']
    d['ref-ephem-tc-2'] = ephem.ix[2,'T']

    inpfn = 'analysis/ttv-lithwick/emcee-lithwick-muprior.py'
    mod = imp.load_source('mod',inpfn)
    d['lithwick-prior-mu1'] = "{:.0f}".format(mod.mu1*1e6)
    d['lithwick-prior-mu1_err'] = "{:.0f}".format(mod.mu1_err*1e6)
    d['lithwick-prior-mu2'] = "{:.0f}".format(mod.mu2*1e6)
    d['lithwick-prior-mu2_err'] = "{:.0f}".format(mod.mu2_err*1e6)

    inpfn = 'analysis/ttv-ttvfast/mnest-ttvfast-muprior.py'
    mod = imp.load_source('mod',inpfn)
    d['ttvfast-tstart'] =  mod.tstart 

    inpfn = 'analysis/ttv-ttvfast/emcee-ttvfast-muprior.py'
    mod = imp.load_source('mod',inpfn)
    d['ttvfast-emcee-nwalkers'] =  mod.nwalkers
    d['ttvfast-emcee-nsteps'] =  mod.nsteps
    d['ttvfast-emcee-nburn'] =  mod.nburn

    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d

    return lines

def val_fit():
    d = OrderedDict()
    chain =  ktwo19.io.load_table('ktwo19_npl=3-ccc',cache=1)

    fmt = OrderedDict()
    fmt['per1'] = "{:.5f}"
    fmt['per2'] = "{:.5f}"
    fmt['per3'] = "{:.0f}"
    fmt['tc1'] = "{:.2f}"
    fmt['tc2'] = "{:.2f}"
    fmt['tc3'] = "{:.0f}"
    fmt['k1'] = "{:.1f}"
    fmt['k2'] = "{:.1f}"
    fmt['k3'] = "{:.1f}"
    fmt['a1'] = "{:.3f}"
    fmt['a2'] = "{:.3f}"
    fmt['a3'] = "{:.2f}"
    fmt['mpsini1'] = "{:.1f}"
    fmt['mpsini2'] = "{:.1f}"
    fmt['mpsini3'] = "{:.0f}"
    fmt['musini1'] = "{:.1f}"
    fmt['musini2'] = "{:.1f}"
    fmt['musini3'] = "{:.1f}"
    fmt['gamma_j'] = '{:.1f}'
    fmt['jit_j'] = '{:.1f}'
    fmt['dvdt'] = '{:.1f}'

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

    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)
    return lines

def val_lithwick():
    """
    Print values associated with the Lithwick fits
    """
    chain = ktwo19.io.load_table('lithwick-emcee-samples',cache=1)
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

def val_ttvfast():
    d = OrderedDict()
    samp = ktwo19.io.load_table('ttvfast-emcee-samples-derived',cache=1)

    fmt = OrderedDict()
    fmt['muppm1'] = "{:.1f}"
    fmt['per1'] = "{:.4f}"
    fmt['secosw1'] = "{:.2f}"
    fmt['sesinw1'] = "{:.2f}"
    fmt['tc1'] = "{:.4f}"
    fmt['muppm2'] = "{:.1f}"
    fmt['per2'] = "{:.4f}"
    fmt['secosw2'] = "{:.2f}"
    fmt['sesinw2'] = "{:.2f}"
    fmt['tc2'] = "{:.4f}"

    # Derived parameters
    fmt['masse1'] = "{:.1f}"
    fmt['masse2'] = "{:.1f}"
    fmt['mr2'] = "{:.2f}"
    fmt['e1'] = "{:.2f}"
    fmt['e2'] = "{:.2f}"
    fmt['prad1'] = "{:.1f}"
    fmt['prad2'] = "{:.1f}"

    fmt['rho1'] = "{:.2f}"
    fmt['rho2'] = "{:.2f}"

    insert_chain_dict(samp, d, fmt) 
    d['e1_p90']  = "{:.2f}".format(samp['e1'].quantile(0.9))
    d['e2_p90']  = "{:.2f}".format(samp['e2'].quantile(0.9))

    # envelope fraction
    fmt = OrderedDict()
    samp = ktwo19.io.load_table('fenv-samples',cache=1)
    fmt['fenv1'] = "{:.0f}"
    fmt['fenv2'] = "{:.0f}"
    fmt['mcore1'] = "{:.1f}"
    fmt['mcore2'] = "{:.1f}"
    fmt['menv1'] = "{:.1f}"
    fmt['menv2'] = "{:.1f}"

    insert_chain_dict(samp, d, fmt) 

    ## Eccentric model

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

def print_table(method):
    if method=="lithwick2-muprior":
        table(mod, samples, fmtd)
    elif method=="lithwick2":
        samples = ktwo19.io.load_samples('lithwick2')
        samples['muppm1'] = samples.mu1 * 1e6
        samples['muppm2'] = samples.mu2 * 1e6
        samples['zmag']= samples.eval('sqrt( rezfree**2 + imzfree**2)') 
        samples['mr2']= samples.eval('mu2 / mu1') 

        import ktwo19.mnest.mnest_lithwick2  as mod
        fmtd = OrderedDict()
        fmtd['per1'] = 5
        fmtd['tc1'] = 4
        fmtd['muppm1'] = 1
        fmtd['per2'] = 4
        fmtd['tc2'] = 4
        fmtd['muppm2'] = 1
        fmtd['ReZfree'] = 3
        fmtd['ImZfree'] = 3
        fmtd['zmag'] = 3
        fmtd['mr2'] = 2
        table(mod, samples, fmtd)
    elif method=='ttvfast-npar=11':
        samples = ktwo19.io.load_samples(method)
        samples['mass2'] = samples['mr2'] * samples['mass1'] 
        samples['masse1'] = samples['mass1'] / M_earth
        samples['masse2'] = samples['mass2'] / M_earth
        samples['e1'] = samples.eval('sqrt(ecosw1**2 + esinw1**2)')
        samples['e2'] = samples.eval('sqrt(ecosw2**2 + esinw2**2)')

        import ktwo19.mnest.ttvfast_npar9  as mod
        fmtd = OrderedDict()
        fmtd['mass1'] = 5
        fmtd['per1'] = 5
        fmtd['ecosw1'] = 3
        fmtd['esinw1'] = 3
        fmtd['inc1'] = 0
        fmtd['node1'] = 0
        fmtd['tc1'] = 4
        fmtd['mr2'] = 2
        fmtd['per2'] = 4
        fmtd['ecosw2'] = 3
        fmtd['esinw2'] = 3
        fmtd['inc2'] = 0
        fmtd['node2'] = 0
        fmtd['tc2'] = 4
        fmtd['stellar_mass'] = 2
        fmtd['masse1'] = 1
        fmtd['masse2'] = 1
        fmtd['e1'] = 3
        fmtd['e2'] = 3
        table(mod, samples, fmtd)
    else:
        raise ValueError, "{} not supported".format(method)
