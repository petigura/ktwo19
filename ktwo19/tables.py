import numpy as np
from ktwo24.config import M_earth
import ktwo24.io
from collections import OrderedDict
def table(mod, samples, format_dict):
    """
    Print quantiles from MCMC fitting.
    """

    # Load in the fixed parameters
    for col, prec in format_dict.iteritems():
        
        if list(samples.columns).count(col)==1:
            samples_to_values(samples, col, prec)
        else:
            
            assert list(mod.parameters_full).count(col)==1, \
                "missing column {}".format(col)
            idx = np.where(mod.parameters_full==col)[0][0]
            val = mod.fixed_val[idx]
            valstr = "{0:.{1:}f}".format(val, prec)
            print r"{{{:s}}}{{ {:s} }}%".format(col, valstr)
         
def samples_to_values(samples, col, prec):
    """
    Convert MCMC samples into a nicely formatted latex string.
    """
    p14, p50, p86 = np.percentile(samples[col],[14,50,86])
    val = p50
    err1 = p86 - p50
    err2 = p14 - p50
    valstr ="{1:.{0:}f}^{{ +{2:.{0:}f} }}_{{ {3:.{0:}f} }}".format(
        prec, val, err1, err2
    )
    print r"{{{:s}}}{{ {:s} }}%".format(col, valstr)
    

def print_table(method):
    if method=="lithwick2-muprior":
        samples = ktwo24.io.load_samples('lithwick2-muprior')
        samples['muppm1'] = samples.mu1 * 1e6
        samples['muppm2'] = samples.mu2 * 1e6
        samples['zmag']= samples.eval('sqrt( ReZfree**2 + ImZfree**2)') 
        samples['mr2']= samples.eval('mu2 / mu1') 

        import ktwo24.mnest.mnest_lithwick2_muprior  as mod
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
    elif method=="lithwick2":
        samples = ktwo24.io.load_samples('lithwick2')
        samples['muppm1'] = samples.mu1 * 1e6
        samples['muppm2'] = samples.mu2 * 1e6
        samples['zmag']= samples.eval('sqrt( ReZfree**2 + ImZfree**2)') 
        samples['mr2']= samples.eval('mu2 / mu1') 

        import ktwo24.mnest.mnest_lithwick2  as mod
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
        samples = ktwo24.io.load_samples(method)
        samples['mass2'] = samples['mr2'] * samples['mass1'] 
        samples['masse1'] = samples['mass1'] / M_earth
        samples['masse2'] = samples['mass2'] / M_earth
        samples['e1'] = samples.eval('sqrt(ecosw1**2 + esinw1**2)')
        samples['e2'] = samples.eval('sqrt(ecosw2**2 + esinw2**2)')

        import ktwo24.mnest.ttvfast_npar9  as mod
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
