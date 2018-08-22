import numpy as np
import radvel
import os
from scipy import optimize

import imp
from ktwo19.config import DATADIR

def posterior():
    default = os.path.join(DATADIR,'../analysis/radvel/','ktwo19_npl=3-default.py')
    P = imp.load_source('test', default)
    for d in dir(P):
        if d[0]!='_':
            exec("{}  = P.{}".format(d,d))


    gphot = ktwo19.io.load_table('phot-gp')
    data = ktwo19.io.load_table('rv',cache=1)
    model = radvel.model.RVModel(params,time_base=time_base)

    # Add in GP parameters
    hnames = [
        'gp_amp', # eta_1; GP variability amplitude
        'gp_explength', # eta_2; GP non-periodic characteristic length
        'gp_per', # eta_3; GP variability period
        'gp_perlength', # eta_4; GP periodic characteristic length
    ]
    params['gp_amp'] = radvel.Parameter(value=5,vary=True)
    params['gp_explength'] = radvel.Parameter(value=gphot.eta2,vary=True)
    params['gp_per'] = radvel.Parameter(value=gphot.eta3,vary=True)
    params['gp_perlength'] = radvel.Parameter(value=gphot.eta4,vary=True)
    like = radvel.likelihood.GPLikelihood(
        model, data.time, data.mnvel,data.errvel, hnames, suffix='_j',
        kernel_name="QuasiPer"
    )

    model = radvel.model.RVModel(params,time_base=2500)
    like.params['gamma_j'] = radvel.Parameter(value=data.mnvel.mean())
    like.params['jit_j'] = radvel.Parameter(value=1)
    like = radvel.likelihood.CompositeLikelihood([like])

    post = radvel.posterior.Posterior(like)    
    post.priors += [radvel.prior.Jeffreys('gp_amp', 0.1, 30.)]
    post.priors += [radvel.prior.Gaussian('gp_explength', gphot.eta2, gphot.eta2_err)]
    post.priors += [radvel.prior.Gaussian('gp_per', gphot.eta3, gphot.eta3_err)]
    post.priors += [radvel.prior.Gaussian('gp_perlength', gphot.eta4, gphot.eta4_err)]
    post.priors += [radvel.prior.Jeffreys('jit_j', 0.01, 10.)]
    post.params['dvdt'].vary=True
    return post

def maximum_a_posteriori():
    post = posterior()
    res = optimize.minimize(
        post.neglogprob_array, post.get_vary_params(), method='Powell',
        options=dict(maxiter=200, maxfev=100000, xatol=1e-8)
    )
    return post


def mcmc():
    post = max_a_posteriori()
    chains = radvel.mcmc(post,nrun=1000,ensembles=3)
    return chains
