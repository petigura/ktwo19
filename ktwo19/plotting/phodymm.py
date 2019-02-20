from chainconsumer import ChainConsumer
import ktwo19.phodymm
import ktwo19.io
from matplotlib.pylab import *

    

def key2tex(k):
    d = {
        'masse2':'$M_{p,b}\, (M_\oplus)$',
        'masse3':'$M_{p,c}\, (M_\oplus)$',
        'ecosw2':'$e_{b}\, \cos \omega_b$',
        'ecosw3':'$e_{c}\, \cos \omega_c$',
        'esinw2':'$e_{b}\, \sin \omega_b$',
        'esinw3':'$e_{c}\, \sin \omega_c$',
    }
    return d[k]

class Plotter(object):
    def __init__(self):
        self.chain = ktwo19.io.load_table('photodyn-samples',cache=1)

    def plot_corner_ecc(self):
        c = ChainConsumer()
        chain = self.chain_without_burnin()
        cols = 'ecosw2 esinw2'.split()
        parameters = [r'$e \cos \omega$',r'$e \sin \omega$']
        c.add_chain(np.array(chain[cols]),parameters=parameters, name="K2-19b")
        cols = 'ecosw3 esinw3'.split()
        c.add_chain(np.array(chain[cols]),parameters=parameters, name="K2-19c")
        c.configure(plot_hists=False)
        c.plotter.plot(figsize=(3,3))
        fig = gcf()
        axL = fig.get_axes()
        setp(axL ,xlim=(-0.301,0.3),ylim=(-0.301,0.3))
        fig.set_tight_layout(True)

    def plot_corner_mass(self):
        c = ChainConsumer()
        chain = self.chain_without_burnin()
        cols = 'masse2 masse3'.split()
        parameters = [key2tex(k) for k in cols]
        c.add_chain(np.array(chain[cols]),parameters=parameters)
        c.configure(plot_hists=False)
        c.plotter.plot(figsize=(3,3))
        fig = gcf()
        axL = fig.get_axes()
        setp(axL ,xlim=(-0.01,40))
        setp(axL,ylim=(-0.01,40))
        fig.set_tight_layout(True)

    def plot_corner_ecosw(self):
        c = ChainConsumer()
        chain = self.chain_without_burnin()
        cols = 'ecosw2 ecosw3'.split()
        parameters = [key2tex(k) for k in cols]
        c.add_chain(np.array(chain[cols]),parameters=parameters)
        c.configure(plot_hists=False)
        c.plotter.plot(figsize=(3,3))
        fig = gcf()
        axL = fig.get_axes()
        setp(axL, xlim=(-0.301,0.3))
        setp(axL, ylim=(-0.301,0.3))
        fig.set_tight_layout(True)

    def plot_corner_esinw(self):
        c = ChainConsumer()
        chain = self.chain_without_burnin()
        cols = 'esinw2 esinw3'.split()
        parameters = [key2tex(k) for k in cols]
        c.add_chain(np.array(chain[cols]),parameters=parameters)
        c.configure(plot_hists=False)
        c.plotter.plot(figsize=(3,3))
        fig = gcf()
        axL = fig.get_axes()
        setp(axL, xlim=(-0.301,0.3))
        setp(axL, ylim=(-0.301,0.3))
        fig.set_tight_layout(True)


    def plot_corner_all(self):
        c = ChainConsumer()
        chain = self.chain_without_burnin()
        cols = 'per1 tc1 inc1 masse1 ror1 '
        cols +='per2 tc2 secosw2 secosw2 inc2 masse2 ror2 ' 
        cols +='per3 tc3 secosw3 secosw3 inc2 Omega3 masse3 ror3 ' 
        cols +='rstar mstar c1 c2'
        cols = cols.split()
        c.add_chain(np.array(chain[cols]),parameters=cols)
        c.configure(plot_hists=False)
        c.plotter.plot(figsize='GROW')
        fig = gcf()
        axL = fig.get_axes()
        fig.set_tight_layout(True)



    def chain_without_burnin(self):
        return self.chain[self.chain.niter > self.nburn /100]

    def plot_chain_trace(self):
        c = ChainConsumer()
        cols = 'per2 tc2 secosw2 sesinw2 per3 tc3 secosw3 sesinw3 Omega3'.split()
        j = 0 
        nchain = 10 
        g = self.chain.groupby('chain')
        for chain, idx in g.groups.iteritems():
            c.add_chain(np.array(self.chain.loc[idx,cols]),parameters=cols)
            j+=1
            if j > nchain:
                break

        c.plotter.plot_walks()

    def plot_corner(self, nburn=10000):
        c = ChainConsumer()
        cols = 'per2 per3 ecosw2 ecosw3 esinw2 esinw3 inc2 inc3 masse2 masse3'.split()
        chain = self.chain_without_burnin()
        c.add_chain(chain,parameters=cols)
        c.configure(sigmas=[0,1,2])
        c.plotter.plot()

    def plot_corner_interesting(self):
        c = ChainConsumer()
        cols = 'ecosw2 ecosw3 esinw2 esinw3 inc2 inc3 Omega3 masse2 masse3'.split()
        chain = self.chain_without_burnin()
        c.add_chain(np.array(chain[cols]),parameters=cols)
        c.configure(sigmas=[0,1,2])
        c.plotter.plot()

def fig_corner(mode):
    fname = 'analysis/photodyn/runs/K2-19_e-uniform_Omega-vary/k2-19.in'
    demcmcfname = 'analysis/photodyn/runs/K2-19_e-uniform_Omega-vary/demcmc_k2-19_massprior.out'
    p = ktwo19.plotting.phodymm.Plotter(fname,demcmcfname)
    p.nburn = 10000
    if mode=='corner-ecc':
        p.plot_corner_ecc()
    elif mode=='corner-mass':
        p.plot_corner_mass()
    elif mode=='corner-ecosw':
        p.plot_corner_ecosw()
    elif mode=='corner-esinw':
        p.plot_corner_esinw()
    elif mode=='corner-all':
        p.plot_corner_all()

import pandas as pd

import matplotlib.gridspec as gridspec

def fig_photodyn_bestfit():
    df = pd.read_csv('analysis/photodyn/runs/K2-19_e-uniform_Omega-vary_2/lc_k2-19_massprior.lcout')
    lcdatafile = "analysis/photodyn/runs/K2-19_e-uniform_Omega-vary/lc_k2-19_massprior.lcout"
    names = 't ig1 flux ig2'.split() 
    fname = 'analysis/photodyn/runs/K2-19_e-uniform_Omega-vary_2/lc_k2-19_massprior.lcout'
    df = pd.read_csv(fname,sep='\s+',names=names)
    lcdata = np.loadtxt(lcdatafile)
    time = lcdata[:,0]
    meas = lcdata[:,1]
    the = lcdata[:,2]
    err = lcdata[:,3]
    nrows = 3
    ncols = 5
    fig,axL = subplots(ncols=nrows, nrows=ncols,figsize=(8,6),sharey=True)
    axL = axL.flatten()
    xc = [
        1980.38, 1984.27, 1988.30, 1996.18, 
        2004.14, 2008.09, 2019.98, 2020.00, 
        2027.90, 2031.91, 2035.85, 2043.74, 
        2043.81, 2051.66, 2055.72]
    buf = 0.25
    i = 0

    row, col = mgrid[0:nrows,0:ncols]
    row = row.flatten()
    col = col.flatten()
    gs0 = gridspec.GridSpec(nrows, ncols)
    for i in range(len(xc)):
        _row = row[i]
        _col = col[i]
        gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[i],height_ratios=[3, 1],hspace=0.0)
    #    gs00.update(vspace=0.05)
        ax1 = plt.subplot(gs00[0])
        ax2 = plt.subplot(gs00[1],sharex=ax1)
        _xc = xc[i]
        sca(ax1)
        errorbar(time-_xc,meas,yerr=err,fmt='.')

        plot(df.t-_xc,df.flux)
        ylim(.9901,1.0019)
        t = ax1.get_xticklabels()
        setp(t, visible=False)
        sca(ax2)
        ylim(-1.5e-3,1.5e-3)

        errorbar(time-_xc,meas - the,yerr=err,fmt='.',ms=1)
        xlim(-buf,+buf)
        if _col!=0:
            t = ax1.get_yticklabels()
            setp(t, visible=False)
            t = ax2.get_yticklabels()
            setp(t, visible=False)
        else:
            setp(ax1,ylabel='Flux')
            setp(ax2,ylabel='Resid')

        if _row!=2:
            t = ax1.get_xticklabels()
            setp(t, visible=False)
            t = ax2.get_xticklabels()
            setp(t, visible=False)
        else:
            setp(ax2,xlabel='Time (days)')
        

        if i>=len(xc):
            setp([ax1,ax2],visible=False)
        
        i+=1

    fig.subplots_adjust(hspace=0.001,wspace=0.001,left=0.1,right=0.97,bottom=0.07,top=0.97)





import numpy as np
import matplotlib.pyplot as plt
import glob
def read_in_tbv(fname):
    data = np.loadtxt(fname)
    n = data[:,0]
    tt = data[:,1]
    sorti = np.argsort(n)
    n = n[sorti]
    tt = tt[sorti]
    return n, tt

def compute_omc(n, tt):
    fit = np.polyfit(n, tt, 1)
    b = fit[0]
    a = fit[1]
    c = b*n + a
    omc = 24.*60.*(tt-c)
    return tt, omc, b

def plot_omc(tt, omc, b, fname=None):
    ax = gca()
    ax.scatter(tt, omc)
    plt.text(0.02, 0.95, "Average Period = %.5lf days" % b, transform=ax.transAxes)
    plt.xlabel('Time (days)')
    plt.ylabel('TTV (minutes)')



    
class TTVPlotter(object):
    """
    
    """


    def __init__(self):
        dirs = glob.glob('K2-19_e-uniform_Omega-vary_mcmc-draws/draws-results/draws/draw*.aei/')
        j = 0 
        self.nsamp = len(dirs)
        self.mod = [[],[]]
        self.obs = [[],[]]
        self.mcolors = ['RoyalBlue','Tomato']

        for i in range(2,4):
            isamp = 0 
            for _dir in dirs:
                f = '{}/tbv00_0{}.out'.format(_dir,i)
                n, tt = read_in_tbv(f)
                if isamp==0:
                    fit = np.polyfit(n, tt, 1)
                    self.b = fit[0]
                    self.a = fit[1]
                    self.c = self.b*n + self.a
                f2 = 'K2-19_e-uniform_Omega-vary_mcmc-draws/ttv_0{}.in'.format(i)
                obs = pd.read_csv(f2,sep='\s+',names=['n','t','terr'])
                obs['omc'] = obs.t - (self.b*obs.n + self.a) 
                obs['omc'] *= 24.
                obs['omc_err'] = obs.terr * 24 

                omc = 24.*(tt-self.c)
                mod = dict(n=n,t=tt,omc=omc)
                mod = pd.DataFrame(mod)
                obs = pd.merge(obs,mod,on='n',suffixes=['_obs','_mod'])
                obs['resid'] = obs.omc_obs - obs.omc_mod
                isamp+=1
                self.mod[j].append(mod)
                self.obs[j].append(obs)

            j+=1

    def plot_ttv(self):
        for j in range(2):
            obs = self.obs[j][0]
            errorbar(obs.t_obs,obs.omc_obs,yerr=obs.omc_err,fmt='o',color='k')
            for isamp in range(self.nsamp):
                mod = self.mod[j][isamp]
                plot(mod.t, mod.omc,color=self.mcolors[j],alpha=0.1)
#            plot(temp.t_mod,resid,color=mcolors[j],alpha=0.1)
            

    def plot_resid(self, j ):
        obs = self.obs[j][0]
        x = np.array(obs.t_mod)
        y = np.zeros_like(obs.t_mod)
        yerr = np.array(obs.omc_err)
        errorbar(x,y,yerr=yerr,fmt='_',color='k',ms=0,zorder=9,mew=2)
        for isamp in range(self.nsamp):
            obs = self.obs[j][isamp]
            plot(obs.t_obs, obs.resid, color=self.mcolors[j],alpha=0.1)

def fig_ttvfit():
    import seaborn as sns
    sns.set_context('paper',font_scale=1.2)
    sns.set_style('ticks')
    gs = gridspec.GridSpec(100, 1)
    fig = figure(figsize=(8,5))
    ax0 = plt.subplot(gs[:60])
    ax1 = plt.subplot(gs[60:80],sharex=ax0)
    ax2 = plt.subplot(gs[80:],sharex=ax0)
    axL = [ax0,ax1,ax2]
    sca(axL[0])
    p = ktwo19.plotting.phodymm.TTVPlotter()
    p.plot_ttv()

    sca(axL[1])
    p.plot_resid(0)
    
    sca(axL[2])
    p.plot_resid(1)
    setp(axL[0],xlim=(1900,5000),ylim = (-12,12))
#    setp([axL[1],axL[2]],ylim = (-0.15,0.15),ylabel='Resid \n (hours)')

    setp([axL[2]],ylim =(-0.15,0.15),ylabel='                  Residuals (hours)')

    setp(axL[2],xlabel='BJD - 2545833')
    setp(axL[0],ylabel='TTV (hours)')

    for ax in axL:
        sca(ax)
        grid()

    for ax in [axL[0],axL[1]]:
        sca(ax)
        t = ax.xaxis.get_ticklabels()
        setp(t,visible=False)

    tight_layout()
    fig.subplots_adjust(hspace=0.001)

#        xlim(1900,4000)
#        ylim(-10,10)





