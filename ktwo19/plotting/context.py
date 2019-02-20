from astropy import constants as c 
import seaborn as sns
import ktwo19.values
sns.set_style('ticks')
sns.set_context('paper',font_scale=1.2)
import subsaturn.lopez
lopi = subsaturn.lopez.LopezInterpolator()
from matplotlib.pylab import *
import pandas as pd
def fig_massradius():

    fig = figure(figsize=(3.0,3.0))
    df = pd.read_csv('planets.csv',comment='#')
    df = df.query('1000 > pl_bmassj / pl_bmassjerr1 > 5 and 100 > pl_radj / pl_radjerr1 > 4 and pl_radjlim==0.0 and pl_bmassjlim==0.0 ')
    df = df[df.pl_discmethod.str.contains('Tr|RV')]
    df = df[~df.pl_name.str.contains('Kepler-128')]
    df['pl_bmasse'] = df['pl_bmassj'] * c.M_jup/ c.M_earth
    df['pl_rade'] = df['pl_radj'] * c.R_jup/ c.R_earth
    df['pl_bmasseerr1'] = df['pl_bmassjerr1'] * c.M_jup/ c.M_earth
    df['pl_radeerr1'] = df['pl_radjerr1'] * c.R_jup/ c.R_earth
    #plot(df.pl_bmasse, df.pl_rade,'.')
    errorbar(df.pl_bmasse,df.pl_rade,xerr=np.array(df.pl_bmasseerr1),yerr=np.array(df.pl_radeerr1),fmt='.',color='gray',elinewidth=1)

    loglog()
    df2 = []
    d = ktwo19.values.val_sys(return_dict=True)
    for i,pl_name in zip([2,3],['b','c']):
        _d = {}
        _d['pl_name'] = "K2-19" + pl_name
        _d['pl_masse'] = float(d['pd-masse%i' % i])
        _d['pl_masseerr'] = float(d['pd-masse%i_err1'%i ])
        _d['pl_rade'] = float(d['pd-prad%i' % i ])
        _d['pl_radeerr'] = float(d['pd-prad%i_err1' % i ])
        df2.append(_d)

    df2 = pd.DataFrame(df2)

    errorbar(df2.pl_masse,df2.pl_rade,xerr=np.array(df2.pl_masseerr),yerr=np.array(df2.pl_radeerr),fmt='.',color='r')
    df2.apply(lambda x : text(x.pl_masse,x.pl_rade,x.pl_name+'  ',weight='bold',ha='right',size='small'),axis=1)

    errorbar(df2.pl_masse,df2.pl_rade,xerr=np.array(df2.pl_masseerr),yerr=np.array(df2.pl_radeerr),fmt='.',color='r')
    xlim(1,1000)
    xlabel('Planet mass (Earth-Masses)')
    ylabel('Planet size (Earth-Radii)')
    ylim(1,30)
    #errorbar(float(d['pd-masse2']),float(d['pd-prad2']),yerr=float(d['pd-prad2_err1']),xerr=float(d['pd-masse2_err1']),fmt='.')
    #errorbar(float(d['pd-masse3']),float(d['pd-prad3']),yerr=float(d['pd-prad3_err1']),xerr=float(d['pd-masse3_err1']),fmt='.')
    tight_layout()

    def plot_massrad(fenv):
        mass = np.logspace(np.log10(3),np.log10(100),1000)
        kw = dict(ls='--',color='blue',lw=0.5) 
        rad = mass2radius(mass,fenv)
        plot(mass,rad,**kw)
        text(mass[-1],rad[-1],"$f_\mathrm{env}$ = %i%%" % (100*fenv),va='center',size='small' )

    plot_massrad(0.1)
    plot_massrad(0.2)
    plot_massrad(0.5)

    xt = [1,3,10,30,100,300,1000]
    xticks(xt,xt)
    yt = [1,2,4,8,16,32]
    xticks(xt,xt)

    gcf().savefig('fig_context-massradius.pdf')

def mass2radius(mass, fenv):
    flux = 70
    time = 5
    comp = 1-fenv
    logcomp = np.log10(comp)
    logflux = np.log10(flux)
    logtime = np.log10(time * 1e9)
    radius = []
    for i in range(len(mass)):
        _mass = mass[i]
        _logmass = np.log10(_mass)
        point_i = [logtime, logflux, logcomp, _logmass]
        logradius = lopi._logradius_interpolate(point_i)
        radius.append(10**logradius)
    radius = np.array(radius)
    return radius


def load_context():

    columns = 'pl_name pl_num pl_radius pl_mass pl_density pl_orbeccen pl_teq pl_fenv st_metfe'.split()
    df = pd.read_csv('data/tab_subsaturns.tex',sep='&',header=None,names=columns,skiprows=1,index_col=None)


    df2 = []
    for i,row in df.iterrows():
        d = {}
        d['pl_name'] = row['pl_name']
        cols = 'pl_radius pl_mass pl_fenv'.split()
        for k in cols:
            s =row[k]
            if s.count('nodata'):
                continue 
            val,err = s.split('^')
            err1, err2 = err.split('_')
            err1 = err1.replace('{','').replace('}','')
            err2 = err2.replace('{','').replace('}','')
            d[k] = val
            d[k+'err1'] = err1
            d[k+'err2'] = err2 
        cols = 'pl_teq st_metfe'.split()
        for k in cols:
            d[k] = row[k]

        df2.append(d)

    df2 = pd.DataFrame(df2)
    df2['st_metfe'] = df2.st_metfe.str.replace('\\','').str.replace('nodata','').str.strip()
    for c in 'pl_radius pl_mass pl_fenv'.split():
        df2[c] = pd.to_numeric(df2[c])
        df2[c+'err1'] = pd.to_numeric(df2[c+'err1'])
        df2[c+'err2'] = pd.to_numeric(df2[c+'err2'])

    for c in 'pl_teq st_metfe'.split():
        df2[c] = pd.to_numeric(df2[c])

    d = ktwo19.values.val_sys(return_dict=True)
    for i,pl_name in zip([2,3],['b','c']):
        _d = {}
        _d['pl_name'] = "K2-19" + pl_name
        _d['pl_mass'] = float(d['pd-masse%i' % i])
        _d['pl_masserr1'] = float(d['pd-masse%i_err1'%i ])
        _d['pl_radius'] = float(d['pd-prad%i' % i ])
        _d['pl_radiuserr1'] = float(d['pd-prad%i_err1' % i ])
        _d['st_metfe'] = 0.06
        _d['pl_fenv'] = float(d['lopez-fenv%i' % i ])
        _d['pl_fenverr1'] = float(d['lopez-fenv%i_err1' % i ])    
        _d['pl_teq'] = float(d['lopez-teq%i' % i ])    

        df2 = df2.append(pd.Series(_d),ignore_index=True)

    df2['pl_name'] = df2['pl_name'].str.strip()
    df2.index = df2['pl_name']
    df2.loc['K2-24 b','pl_mass'] = 19
    df2.loc['K2-24 b','pl_masserr1'] = 2
    df2.loc['K2-24 b','pl_masserr2'] = -2
    df2.loc['K2-24 b','pl_fenv'] = 26
    df2.loc['K2-24 b','pl_fenverr1'] = 3
    df2.loc['K2-24 b','pl_fenverr2'] = -3
    df2.loc['K2-24 c','pl_mass'] =  15 
    df2.loc['K2-24 c','pl_masserr1'] =  2 
    df2.loc['K2-24 c','pl_masserr2'] =  -2 
    df2.loc['K2-24 c','pl_fenv'] = 52
    df2.loc['K2-24 c','pl_fenverr1'] = 4
    df2.loc['K2-24 c','pl_fenverr2'] = -4

    # Adding in brady planet 
    _d = {}
    _d['pl_name'] = "Kepler-1656b"
    _d['pl_mass'] = 48.6
    _d['pl_masserr1'] = 3.8
    _d['pl_radius'] = 5.02
    _d['pl_radiuserr1'] = 0.5
    _d['st_metfe'] = 0.19 
    df2 = df2.append(pd.Series(_d),ignore_index=True)
    return df2

def fig_smetmass():
    df2 = load_context()
    figure(figsize=(3.0,3.0))
    sns.set_style('ticks')
    sns.set_context('paper',font_scale=1.2)
    df2['pl_mcore'] = df2.eval('pl_mass - 0.01 *pl_fenv * pl_mass')
    semilogy()
    
    xk = 'st_metfe'
    yk = 'pl_mass'
    yerrk ='pl_masserr1'
    #yk = 'pl_mcore'

    kw = dict(fmt='o',ms=5)
    errorbar(df2[xk],df2[yk],yerr=df2[yerrk],color='gray',**kw)
    df3 = df2[df2.pl_name.str.contains('K2-19')]
    errorbar(df3[xk],df3[yk],yerr=df3[yerrk],color='red',**kw)
    df3.apply(
        lambda x : text(x[xk],x[yk],x.pl_name+'  ',weight='bold',
                        ha='right',size='small'),axis=1
    )
    xlim(-0.3,0.5) 
    ylim(3,100) 
    xlabel('[Fe/H] (dex)')
    ylabel('Planet Mass (Earth-masses)')
    yt = [1,3,10,30,100]
    yticks(yt,yt)
    tight_layout()

def fig_teqfenv():
    df2 = load_context()
    figure(figsize=(3.0,3.0))
    sns.set_style('ticks')
    sns.set_context('paper',font_scale=1.2)

    xk = 'pl_teq'
    yk = 'pl_fenv'
    yerrk ='pl_fenverr1'
    kw = dict(fmt='o',ms=5)
    errorbar(df2[xk],df2[yk],yerr=df2[yerrk],color='gray',**kw)
    df3 = df2[df2.pl_name.str.contains('K2-19')]
    errorbar(df3[xk],df3[yk],yerr=df3[yerrk],color='red',**kw)
    df3.apply(
        lambda x : text(x[xk],x[yk],x.pl_name,weight='bold',
                        ha='left',size='small'),axis=1
    )

    df3 = df2[df2.pl_name.str.contains('K2-24 c')]
    df3.apply(
        lambda x : text(x[xk],x[yk],x.pl_name,weight='bold',
                        ha='left',size='small'),axis=1
    )

    xlabel('Equilibrium Temp (K)') 
    ylabel('Envelope Fraction (%)')
    xlim(250,1750) 
    ylim(0,100) 
    tight_layout()
