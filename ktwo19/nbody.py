import rebound
import numpy as np
from astropy import constants as c
import radvel.orbit
import ktwo19.io
from matplotlib import pylab as plt
import pandas as pd

def setupSimulation(row):
    sim = rebound.Simulation()
    sim.t = 2020.0 
    sim.G = 2.9591220363e-4 # in Msun-AU-day units
    sim.integrator = "ias15" 
    sim.add(m=row.mstar) # K2-19 
    
    # Add K2-19b
    m = row.masse2 * c.M_earth / c.M_sun # Solar masses
    P = row.per2
    e = row.e2
    omega = row.omega2
    Omega = row.Omegarad2
    inc = row.incrad2
    T = row.tp2 # Time of periastron
    sim.add(m=m, P=P, e=e, omega=omega, Omega=Omega, inc=inc, T=T)

    # Add K2-19c
    m = row.masse3 * c.M_earth / c.M_sun # Solar masses
    P = row.per3
    e = row.e3
    omega = row.omega3
    Omega = row.Omegarad3
    inc = row.incrad3
    T = row.tp3 # Time of periastron
    sim.add(m=m, P=P, e=e, omega=omega, Omega=Omega, inc=inc, T=T)
    sim.move_to_com()

    return sim

def run_simulation(irow):
    p = ktwo19.plotting.phodymm.Plotter()
    p.nburn=10000
    df = p.chain_without_burnin()
    
    # map irow to random irow
    row = df.sample(1, random_state=irow)
    row = row.iloc[0]

    sim = setupSimulation(row)
    sim.status()
    sim.exit_min_distance = 0.01
    archivefn = "nbody/archive-{}.bin".format(irow)
    
    sim.automateSimulationArchive(
        archivefn,interval=10,deletefile=True
    )
    tstart = sim.t
    tstop = 50000
    Noutputs = 1000

    times = np.linspace(tstart, tstop, Noutputs)
    a1 = np.zeros(Noutputs)
    e1 = np.zeros(Noutputs)
    M1 = np.zeros(Noutputs)
    n1 = np.zeros(Noutputs)
    w1 = np.zeros(Noutputs)
    pomega1 = np.zeros(Noutputs)
    a2 = np.zeros(Noutputs)
    e2 = np.zeros(Noutputs)
    M2 = np.zeros(Noutputs)
    n2 = np.zeros(Noutputs)
    w2 = np.zeros(Noutputs)
    pomega2 = np.zeros(Noutputs)

    distances = np.zeros(Noutputs)
    ps = sim.particles # ps is now an array of pointers. It will update as the simulation runs.
    t_print = tstart
    t_since_print = 0
    t_print_step = 100000
    try:
        for i, time in enumerate(times):
            sim.integrate(time, exact_finish_time=0)
            a1[i] = ps[1].a
            e1[i] = ps[1].e
            M1[i] = ps[1].M  # mean anomaly
            n1[i] = ps[1].n  # mean motion
            w1[i] = ps[1].omega # argument 
#            pomega1[i] = ps[1].pomega # longitude
            pomega1[i] = w1[i] + ps[1].Omega
            a2[i] = ps[2].a
            e2[i] = ps[2].e
            M2[i] = ps[2].M  # mean anomaly
            n2[i] = ps[2].n  # mean motion
            w2[i] = ps[2].omega # argument
#            pomega2[i] = ps[2].pomega # longitude
            pomega2[i] = w2[i] + ps[1].Omega
            
            if t_since_print > t_print_step:
                t_print = time
                t_since_print = 0
                print time
            else:
                t_since_print = time - t_print

            dp = ps[1] - ps[2]   # Calculates the coponentwise difference between particles
            distances[i] = np.sqrt(dp.x*dp.x + dp.y*dp.y + dp.z*dp.z)
    except rebound.Encounter as error:
        print(error)

    # confused... should I be plotting delta pomega? It seems like
    # there's some sign flipping going on with poemga
    dw = w2 - w1
    dw = np.arctan2(np.sin(dw), np.cos(dw))

    import seaborn as sns
    sns.set_context('paper',font_scale=0.9)
    fig1 = plt.figure(figsize=(8,2.7))
    ax = plt.subplot(131)
    ax.set_xlabel("BJD - 2454833")
    ax.set_ylabel("Eccentricity")
    plt.plot(times, e1, label='K2-19b',color='RoyalBlue');
    plt.plot(times, e2, label='K2-19c',color='Tomato');
    plt.legend()
    plt.ylim(0,0.3)

    ax = plt.subplot(132)
    ax.set_xlabel(r"$\Delta \varpi$ (deg)")
    ax.set_ylabel("Eccentricity")
    plt.plot(np.rad2deg(dw), e1,'.',ms=2,label='K2-19b',color='RoyalBlue');
    plt.plot(np.rad2deg(dw), e2,'.',ms=2,label='K2-19c',color='Tomato');
    plt.xlim(-180,180)
    plt.ylim(0,0.3)
    plt.legend()

    ax = plt.subplot(133)
    ax.set_xlabel("BJD - 2454833")
    ax.set_ylabel("M1 - M2 [AU]")
    x = 3 * M2
    y = 2 * M1

    #Lambda = np.arctan2(np.sin(x-y), np.cos(x-y))
    plt.plot(times,3 * n2 - 2 * n1, label='$(3 n_c - 2 n_b)$')
    plt.ylabel(r'Angular Velocity (rad/day)')
    #(black), $\dot{\omega}_1$ (rad/day) (red), ')
    plt.tight_layout(True)

    #omegadot1 = (w1[1:] - w1[:-1]) / (times[1:] - times[:-1])
    #omegadot2 = (w2[1:] - w2[:-1]) / (times[1:] - times[:-1])
    omegadot1 = (pomega1[1:] - pomega1[:-1]) / (times[1:] - times[:-1])
    omegadot2 = (pomega2[1:] - pomega2[:-1]) / (times[1:] - times[:-1])

    #Delta = n1/n2 * 2/3 - 1
    plt.plot(times[0:-1],omegadot1,label=r'$\dot{\varpi}_b$')
    plt.plot(times[0:-1],omegadot2,label=r'$\dot{\varpi}_c$')
    plt.ylim(-0.01,0.01)
    plt.legend()

    fig2, axL = plt.subplots(ncols=2,figsize=(9,4),sharey=True)
    ax = axL[0]
    sa = rebound.SimulationArchive(archivefn)
    sim = sa.getSimulation(tstart)
    rebound.plotting.OrbitPlotOneSlice(sim, ax, axes="xz",periastron=True, trails=True )
    plt.grid()
    
    ax = axL[1]
    sim = sa.getSimulation(tstop)
    rebound.plotting.OrbitPlotOneSlice(sim, ax, axes="xz",periastron=True, trails=True )
    plt.grid()

    plt.setp(
        axL, xlim=(-0.15,0.15), ylim=(-0.15,0.15), 
        xlabel='x [AU]',ylabel='z [AU]'
    )
    fig1.savefig('nbody/fig_nbody-elements-{}.pdf'.format(irow))
    fig2.savefig('nbody/fig_nbody-orbits-{}.pdf'.format(irow))




