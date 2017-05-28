from matplotlib.pylab import *
import ttv.plotting
import pandas as pd
import ktwo19.io
import numpy as np
from ktwo19.config import AU_per_day, M_earth
import re
import corner
import copy

ephem = ktwo19.io.load_ephem().copy()
times = ktwo19.io.load_times()
times = times.dropna(subset=['inst'])
times = times[times.inst.str.contains('K2|FLWO|TRAPPIST|MuSCAT|Spitzer')]
times = times.sort_values(by='tc')
times.index = [times.i_planet,times.i_epoch]
levels = [1 - exp(-0.5*x**2) for x in range(1,3)] # levels for corner plots

def sample_resid(ttvmod, chains, i_planet, **kwargs):
    nchain = chains.shape[0]

    t = times.ix[i_planet]
    per = ephem.ix[i_planet,'per']
    T = ephem.ix[i_planet,'T']
    i_epoch = t.i_epoch
    for i_chain in range(nchain):
        res = ttvmod.compute(chains[i_chain])
        positions, rv = ttv.ttvfastmod.parse_ttv_results(res, ttvmod.nplanet)
        times_mod = pd.DataFrame(positions)
        times_mod.index = [times_mod.i_planet,times_mod.i_epoch]    

        tmod = times_mod.ix[i_planet]
        per = ephem.ix[i_planet,'per']
        T = ephem.ix[i_planet,'T']
        resid = (t.tc - tmod.tc).dropna()
        errorbar(np.array(t.tc), resid, yerr=t.tc_err,**kwargs)


def sample_dt(ttvmod, chains, i_planet, **kwargs):
    """Calculate times between transits
    """
    nchain = chains.shape[0]

    for i_chain in range(nchain):
        res = ttvmod.compute(chains[i_chain])
        positions, rv = ttv.ttvfastmod.parse_ttv_results(res, ttvmod.nplanet)
        times_mod = pd.DataFrame(positions)
        times_mod.index = [times_mod.i_planet,times_mod.i_epoch]    

        tmod = times_mod.ix[i_planet]
        time = np.array(tmod.tc)
        dt = time[1:] - time[:-1]
        plot(time[:-1], dt, **kwargs)

def sample_rv(ttvmod, chains, rv_times, **kwargs):
    """Sample the array and compute stellar rvs 
    """
    nchain = chains.shape[0]

    assert min(rv_times) >= ttvmod.tstart +1, "rvtimes out of bounds"
    assert max(rv_times) <= ttvmod.tstop -1, "rvtimes out of bounds"
    for i_chain in range(nchain):
        res = ttvmod.compute(chains[i_chain], rv_times=rv_times)
        times_mod, rv_mod = ttv.ttvfastmod.parse_ttv_results(
            res, ttvmod.nplanet 
            )

        rv_mod *= AU_per_day
        plot(rv_times,rv_mod, **kwargs)

def errorbar_omc(i_planet, **kwargs):
    t = times.ix[i_planet]
    per = ephem.ix[i_planet,'per']
    T = ephem.ix[i_planet,'T']
    ttv.plotting.errorbar_omc(t.i_epoch, t.tc, per, T, **kwargs)

def plot_ttv_samples(res, samples, nsamp=50):

    ttvmod = res.ttvmod
    ttvmod.tstop = 10000
    samp = np.array(samples.sample(nsamp, replace=False))
    def plot_ttvs():
        kw1 = dict(color='red',alpha=0.1)
        kw2 = dict(color='blue',alpha=0.1)

        ttv.plotting.plot_omc(i_epoch1, tc1.T, per, T , **kw1)
        per, T = ephem.ix[2]
        ttv.plotting.plot_omc(i_epoch2, tc2.T, per, T , **kw2)

        sample_omc( ttvmod, samp, 1,**kw1)
        sample_omc( ttvmod, samp, 2,**kw2)
        errorbar_omc(1,fmt='or')
        errorbar_omc(2,fmt='ob')
        ylabel('TTV (days)')
        xlabel('BJD - 2454833')
        legend()

    
    plot_ttvs()
    xlim(2000,3000)
    ylim(-.2,.2)
    figure(figsize=(6,2))

    figure(figsize=(6,2))
    plot_ttvs()
    xlim(2000,6000)
    ylim(-1,1)


def plot_ttv_bestfit(res):
    ttvmod = res.ttvmod
    def plot_ttvs():
        kw1 = dict(color='red',linestyle='--',label='b guess')
        kw2 = dict(color='blue',linestyle='--',label='c guess')
        sample_omc( ttvmod, res.p_restrict0[np.newaxis,:], 1,**kw1)
        sample_omc( ttvmod, res.p_restrict0[np.newaxis,:], 2,**kw2)
        errorbar_omc(1,fmt='or')
        errorbar_omc(2,fmt='ob')
        kw1 = dict(alpha=1,color='red',linestyle='-',label='b fit')
        kw2 = dict(alpha=1,color='blue',linestyle='-',label='c fitguess')
        sample_omc( ttvmod, res.p_restrict1[np.newaxis,:], 1,**kw1)
        sample_omc( ttvmod, res.p_restrict1[np.newaxis,:], 2,**kw2)
        legend()

    plot_ttvs()
    xlim(2000,3000)
    ylim(-.2,.2)

    figure(figsize=(6,2))
    kw1 = dict(fmt='o',color='red',label='b resid')
    kw2 = dict(fmt='o',color='blue',label='c resid')

    sample_resid( ttvmod, res.p_restrict1[np.newaxis,:], 1,**kw1)
    sample_resid( ttvmod, res.p_restrict1[np.newaxis,:], 2,**kw2)
    legend()
    
    figure(figsize=(6,2))
    plot_ttvs()
    xlim(2000,6000)
    ylim(-1,1)

def transit_insets():
    """
    Plot transit fits to the photometry
    """

    s = """
    i_planet itrans
    1 0
    2 0
    1 1
    1 2
    2 1
    1 3
    """
    from cStringIO import StringIO as sio
    from ttv.plotting import transit_insets
    from k2phot.config import bjd0
    df = pd.read_csv(sio(s),sep=' ')
    transL = []

    for i, row in df.iterrows():
        picklefn = 'transit-times/transit_pipe=paper_starname=K2-24_i-planet={i_planet}_itrans={itrans:04d}.pkl'.format(**row)
        print picklefn
        trans = pickle.load(open(picklefn))
        transL += [trans]

    clf()
    lc = pd.read_csv('transit-times/photometry_pipe=paper_starname=K2-24.csv')

    figure(figsize=(7.5,3.5))
    transit_insets(lc,transL)

    sca(gcf().get_axes()[0])
    plt.xlabel('BJD - %i' % bjd0)
    plt.ylabel('Normalized Flux')
    gcf().savefig('paper/fig_k2phot_transits.pdf')

def ttv_rv_rvfits():
    textent = 10 * 365 # plot 10 years of models
    xlim(2000,3000)
    ylim(-20,20)


#from ttv.plotting import zoom_effect01
from matplotlib.ticker import MaxNLocator
npsamp = 100

def rv_insets():
    import ktwo24.mnest.emcee_ttvfast_npar11 as mod
    kw = dict(alpha=0.1,color='r')
    prop = dict(fill=None, linestyle=':',edgecolor='k',lw=0.7)

    samples = pd.read_hdf(
        'analysis/emcee-ttvfast-npar=11_nburn=1000/chain.hdf','chain'
    )

    ttvmod = mod.ttvmod
    samples = samples[(samples.mass1 / M_earth) < 30]
    samp = np.array(samples.sample(npsamp))
    ttvmod.tstart = 2000
    ttvmod.tstop = ttvmod.tstart + 10 * 365

    rv_timesi = linspace(mod.ttvmod.tstart+1,mod.ttvmod.tstop-1,10000)

    nrows = 2
    ncols = 2
    fig = figure(figsize=(7.5,6))
    ax_full = subplot2grid((nrows,ncols), (0,0), colspan=ncols, rowspan=nrows-1)
    ax_insets = [subplot2grid((nrows,ncols), (1,col),) for col in range(ncols)]

    def _plot():
        plot(mod.rv_times, mod.rv,'.')
        sample_rv(mod.ttvmod, samp, list(rv_timesi),**kw)

    width = 140
    plt.sca(ax_full)
    _plot()
    xlim(2200,3000)

    def plot_inset(ax, xstart):
        sca(ax)
        _plot()
        xl = xstart,xstart+width
        xlim(*xl)
        zoom_effect01(ax, ax_full,  xl[0], xl[1], **prop)
        ax.xaxis.set_major_locator(MaxNLocator(4))

    plot_inset(ax_insets[0], 2335)
    grid()
    plot_inset(ax_insets[1], 2700)
    grid()

#    plt.setp(ax_insets,ylim=(0.995,1.001))
#    [plt.setp(ax.xaxis,visible=False) for ax in ax_insets]
    [plt.setp(ax.yaxis.get_ticklabels(),visible=False) for ax in ax_insets[1:]]
    plt.gcf().set_tight_layout(True)
    setp([ax_full,ax_insets[0]],ylabel='RV (m/s)')
    setp([ax_full,ax_insets[0]],xlabel='BJD - 2454833')



ebarkw = dict(mec='w',mew=1,ms=5)

def ttv_samples_ttvfast():
    i_epoch1, tc1 = ktwo19.io.model_from_samples('ttvfast-npar9',1)
    i_epoch2, tc2 = ktwo19.io.model_from_samples('ttvfast-npar9',2)
    rectkw = dict(color='none',ec='blue')
    def plot_ttvs():
        kw1 = dict(color='red',alpha=0.1)
        kw2 = dict(color='blue',alpha=0.1)

        per, T = ephem.ix[1]
        ttv.plotting.plot_omc(i_epoch1, tc1.T, per, T , **kw1)
        per, T = ephem.ix[2]
        ttv.plotting.plot_omc(i_epoch2, tc2.T, per, T , **kw2)
        times.index = times.i_planet
        ktwo19.plotting.errorbar_omc(
            1,yerr=times.ix[1,'tc_err']*5, fmt='or',label='K2-19b TTVs',
            **ebarkw
        )
        ktwo19.plotting.errorbar_omc(
            2,yerr=times.ix[2,'tc_err']*5, fmt='ob',label='K2-19c TTVs',
            **ebarkw)
    
    fig, axL = _times_insets_provision_figure()
    for ax in axL:
        sca(ax)
        plot_ttvs()
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))

        
    sca(axL[0])
    x1 = 1900
    x2 = 3650
    xlim(x1,x1+x2)
    ylim(-1,1)
    legend()
    times_insets()
    setp([axL[0],axL[1],axL[2],axL[5]],ylabel='TTV (days)')
    setp([axL[0]],xlabel='BJD-2454833')
    fig.set_tight_layout(True)

def _times_insets_provision_figure():
    fig = figure(figsize=(6,6))
    rc('font',size=8)
    ncols = 3
    nrows = 4
    ax1 = subplot2grid((nrows,ncols), (0,0), colspan=ncols, rowspan=1)
    ax2 = subplot2grid((nrows,ncols), (1,0), colspan=ncols, rowspan=1)
    ax3 = subplot2grid((nrows,ncols), (2,0),)
    ax4 = subplot2grid((nrows,ncols), (2,1),)
    ax5 = subplot2grid((nrows,ncols), (2,2),)
    ax6 = subplot2grid((nrows,ncols), (3,0),)
    ax7 = subplot2grid((nrows,ncols), (3,1),)

    axL = fig.get_axes()
    return fig, axL 


figlabelprop = prop=dict(weight='bold',size='large')
def times_insets(rectkw={}):
    _times = times.copy()

    rectkw = dict(color='none',ec='blue')
    axL = gcf().get_axes()

    inset(axL[0], axL[1], (2500,0), 1500, 0.6, label='b', rectkw=rectkw)

    at = AnchoredText('a', frameon=True, loc=2, prop=figlabelprop)
    axL[0].add_artist(at)

    xy = _times.iloc[:6].mean()
    xy = xy.tc, 0

    #import pdb;pdb.set_trace()
    ephem='linear'
    if ephem=='linear':
        w=140
        h=0.08
        
        xy = (mean(_times.query('inst=="K2"').tc), -0.01)
        inset(axL[1], axL[2], xy, w, h, label='c',rectkw=rectkw) # transit

        temp = _times[_times.inst.str.contains('FLWO|TRAPPIST|MuSCAT')]
        xy = (mean(temp.tc), -0.00)
        inset(axL[1], axL[3], xy, w, h, label='e',loc='below',rectkw=rectkw) 

        temp = _times[_times.inst.str.contains('Spitzer')]
        xy = (mean(temp.tc), +0.03)
        inset(axL[1], axL[4], xy, w, h, label='d', rectkw=rectkw) # transit



def inset(ax1,ax2, center, width, height, label=None, loc='above', rectkw={}):
    """Draw an inset

    Draws an rectangle on ax1 with a given center, width and height. Sets the limits of ax2 to correspond to ax1 rectangle

    Args:
        ax1: big axis
        ax2: detail axis
    """
    xy = (center[0]-0.5*width, center[1]-0.5*height)
    sca(ax1)
    p = matplotlib.patches.Rectangle(xy, width, height, **rectkw)
    ax1.add_artist(p)
    sca(ax2)
    xlim(center[0]-0.5*width, center[0]+0.5*width)
    ylim(center[1]-0.5*height, center[1]+0.5*height)

    if label is not None:
        if loc=="above":
            x = center[0]
            y = center[1] + 0.5 * height
            va = 'bottom'
        elif loc=="below":
            x = center[0]
            y = center[1] - 0.5 * height
            va = 'top'
        
        sca(ax1)
        text(x,y,label, fontweight='bold', va=va)

        sca(ax2)
        at = AnchoredText(label, frameon=True, loc=2, prop=figlabelprop )
        ax2.add_artist(at)

from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

from matplotlib.transforms import Bbox, TransformedBbox, \
    blended_transform_factory

from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
    BboxConnectorPatch


def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = prop_lines.copy()
        prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    p = BboxConnectorPatch(bbox1, bbox2,
                           # loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           **prop_patches)
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect01(ax1, ax2, xmin, xmax, **kwargs):
    """
    ax1 : the main axes
    ax1 : the zoomed axes
    (xmin,xmax) : the limits of the colored area in both plot axes.

    connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
    be marked.  The keywords parameters will be used ti create
    patches.
    """

    trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
    trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

    bbox = Bbox.from_extents(xmin, 0, xmax, 1)

    mybbox1 = TransformedBbox(bbox, trans1)
    mybbox2 = TransformedBbox(bbox, trans2)

    prop_patches = kwargs.copy()
#    prop_patches["ec"] = "none"
#    prop_patches["alpha"] = 0.1

    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=2, loc2a=3, loc1b=1, loc2b=4,
                     prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)
    return c1, c2, bbox_patch1, bbox_patch2, p


texstr = {
    'stellar_mass':'$M_\star$',
    'mass1':'$M_{P,1}$',
    'per1':'$P_1$',
    'esinw1':'$e\, \sin\, \omega_1$',
    'ecosw1':'$e\, \cos\, \omega_2$',
    'tc1':'$T_{c,1}$',
    'mu1':'$\mu_1$',
    'mass2':'$M_{P,2}$',
    'per2':'$P_2$',
    'esinw2':'$e\, \sin\, \omega_2$',
    'ecosw2':'$e\, \cos\, \omega_2$',
    'tc2':'$T_{c,2}$',
    'mu2':'$\mu_2$',
    'ReZfree':'$\mathrm{Re}(Z_\mathrm{free})$',
    'ImZfree':'$\mathrm{Im}(Z_\mathrm{free})$',
    'mr2':'$M_{P,2}  / M_{P,1}$'
}

format_dict = dict(
    stellar_mass=2, 
    per=5,
    mass=1,
    ecosw=2, 
    esinw=2, 
    inc=0, 
    arg=0, 
    node=0, 
    tc=5,
    mr=2,
    mu=6,
    ReZfree=3,
    ImZfree=3,
)

def corner_ttv_rv_eleven():
    npsamp = 1e5
    samples = ktwo24.io.load_samples('ttvfast-npar=11')
    samples = samples.sample(npsamp)
    cols = 'per1 tc1 per2 tc2 mass1 ecosw1 esinw1 mr2 ecosw2 esinw2 stellar_mass'.split()
    corner_plot(samples,cols,plot_datapoints=False)

def corner_ttv_rv_six():
    npsamp = 1e5
    samples = ktwo24.io.load_samples('ttvfast-npar=11')
    samples = samples.sample(npsamp)
    samples['mass2'] = samples['mass1'] * samples['mr2']
    samples['mass1'] /= M_earth
    samples['mass2'] /= M_earth
    samples = samples.query('mass1 < 30')
    samples = samples.sample(npsamp)
    cols = 'mass1 mass2 ecosw1 esinw1 ecosw2 esinw2'.split()
    corner_plot(samples,cols,plot_datapoints=False)

import corner
import pymultinest


def corner_lithwick_mu_three():
    """
    Corner plot showing the lithwick parameters
    """
    df = ktwo24.io.load_samples('lithwick2-muprior')
    cols = 'mu1 mu2 ReZfree ImZfree'.split()
    corner_plot(df[cols],cols,bins=50,smooth=True,plot_datapoints=False)

def corner_lithwick_mu_eight():
    """
    Corner plot showing the lithwick parameters
    """
    df = ktwo24.io.load_samples('lithwick2-muprior')
    cols = 'per1 tc1 per2 tc2 mu1 mu2 ReZfree ImZfree'.split()
    corner_plot(df[cols],cols,bins=50,smooth=True,plot_datapoints=False)


def corner_lithwick_mu_three2():
    """
    Corner plot showing the lithwick parameters
    """
    df = ktwo24.io.load_samples('lithwick2')
    cols = 'mu1 mu2 ReZfree ImZfree'.split()
    corner_plot(df[cols],cols,bins=100,plot_datapoints=False)

def corner_lithwick_mu_eight2():
    """
    Corner plot showing the lithwick parameters
    """
    df = ktwo24.io.load_samples('lithwick2')
    cols = 'per1 tc1 per2 tc2 mu1 mu2 ReZfree ImZfree'.split()
    corner_plot(df[cols],cols,bins=100,plot_datapoints=False)




def corner_lithwick_three():
    """
    Corner plot showing the lithwick parameters
    """

    def _set_limits(axL):
        setp(axL[0,0],xlim=(-1,100))
        setp(axL[1,0],xlim=(-1,100),ylim=(-1,100))
        setp(axL[2,0],xlim=(-1,100),ylim=(0,0.15))
        setp(axL[3,0],xlim=(-1,100),ylim=(0,0.15))

        setp(axL[1,1],xlim=(-1,100),)
        setp(axL[2,1],xlim=(-1,100),ylim=(0,0.15))
        setp(axL[3,1],xlim=(-1,100),ylim=(0,0.15))

        setp(axL[2,2],xlim=(0.0,0.15))
        setp(axL[3,2],xlim=(0.0,0.15),ylim=(0,0.15))
        setp(axL[3,2],xlim=(0.0,0.15))

    df = ktwo24.io.load_samples('lithwick')
    df = df.query('ReZfree < 0.15') 
    df = df.query('ImZfree < 0.15') 
    cols = 'mu1 mu2 ReZfree ImZfree'.split()
    corner_plot(df[cols],cols,bins=100,hist_kwargs=dict(normed=True))

    fig = gcf()
    axL = np.array(fig.get_axes()).reshape(len(cols),-1)
    _set_limits(axL)

    df = ktwo24.io.load_samples('lithwick-mu-prior')
    cols = 'mu1 mu2 ReZfree ImZfree'.split()
    corner_plot(df[cols],cols,bins=10,fig=fig,hist_kwargs=dict(normed=False),contourf_kwargs=dict(cmap=cm.hot,colors=None,vmin=1e-1))


    axL = np.array(fig.get_axes()).reshape(len(cols),-1)
    _set_limits(axL)

def corner_plot(samples,cols,**kwargs):
    _texstr  = copy.copy(texstr)
    titles = []
    for i, col in enumerate(cols):
        p = col
        _prec = format_func(p)
        fmt = format_func(p)
        if list(cols).count(p) > 0:
            p14, p50, p86 = np.percentile(samples[p],[14,50,86])
            val = p50
            err1 = p86 - p50
            err2 = p14 - p50
            valstr ="{1:.{0:}f}^{{ +{2:.{0:}f} }}_{{ {3:.{0:}f} }}".format(
                _prec, val, err1, err2
            )

        titles.append("{} = ${}$".format(texstr[col],valstr))

    for k in 'per1 per2 tc1 tc2'.split():
        if cols.count(k)==0:
            continue
        smed = "{:.4f}".format(samples[k].median())
        med = float(smed)
        samples[k] -= med
        _texstr[k] = "${} - {}$".format(_texstr[k][1:-1], smed)

    for k in 'mu1 mu2'.split():
        if cols.count(k)==0:
            continue
        smed = "{:.4f}".format(samples[k].median())
        med = float(smed)
        samples[k] *= 1e6
        _texstr[k] = r"${} \times 10^6$".format(_texstr[k][1:-1], smed)


    labels = [_texstr[col] for col in cols]

    rc('font',size=14)
    fig = corner.corner(
        samples[cols],labels=labels,levels=levels,  
        fill_contours=True,show_titles=False,use_math_text=True,
        **kwargs
    )
    axL = fig.get_axes()
    axL = np.array(axL).reshape(len(cols),-1)
    for i, col in enumerate(cols):
        axL[i,i].set_title(titles[i],loc='left',size='large')

    for ax in axL.flatten():
        setp(ax.get_xticklabels(),rotation=35)
        setp(ax.get_yticklabels(),rotation=35)


def format_func(p):
    p = str(p)
    for k in format_dict.keys():
        nmatches = len(re.findall('^'+k, p))
        if nmatches==0:
            continue
        elif nmatches==1:
            return format_dict[k]
        else:
            assert False, "problem"
    
    # If algo makes it here, there is a problem
    assert False, "missing format key"


def dumpbell():
    samples = ktwo24.io.load_samples('lithwick-bootstrap')
    mr2 = samples.eval('-ImV1/ImV2*per2/per1*(1.0/2.0)**(1.0/3.0)')
    samples['mr2'] = mr2

    desc = samples.describe()
    samples = samples.sample(100,random_state=0)
    print "mr2 = {} \pm {}".format(mr2.mean(),mr2.std())
    figure(figsize=(3.5,3.5))
    plot(samples.ReV1,samples.ImV1,'.',label='K2-24b')
    plot(samples.ReV2,samples.ImV2,'.',label='K2-24c')
    legend()
    x = samples['ReV1 ReV2'.split()].median()
    y = samples['ImV1 ImV2'.split()].median()
    plot(x,y)
    axhline(0,color='k',linestyle='--')
    axvline(0,color='k',linestyle='--')
    plt.axis('equal')
    xlim(-0.25,0.25)
    ylim(-0.25,0.25)
    xlabel('$\mathrm{Re}(V)$ (days)')
    ylabel('$\mathrm{Im}(V)$ (days)')
    #figure()
    #hist(mr2)
    
