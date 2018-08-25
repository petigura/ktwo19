import rebound
import numpy as np
from astropy import constants as c
import radvel.orbit
import ktwo19.io
from matplotlib import pylab as plt

def setupSimulation(i):
    df = ktwo19.io.load_table('mills-photodyn-samples',cache=1)
    row = df.iloc[i]
    sim = rebound.Simulation()
    sim.t = 2020.0 
    sim.G = 2.9591220363e-4 # in AU-day units
    sim.integrator = "ias15" # IAS15 is the default integrator, so we don't need this line

    Mstar = 0.88
    sim.add(m=0.88) # K2-19 
    
    # Add K2-19b
    m = row.masse2 * c.M_earth / c.M_sun # Solar masses
    per = row.per2
    a = (per**2 * sim.G * (Mstar + m) / 4 / np.pi**2)**(1/3.)
    e = row.e2
    omega = row.w2
    Omega = row.Omegarad2
    inc = row.incrad2
    T = row.tp2 # Time of periastron
    sim.add(m=m, a=a, e=e, omega=omega, Omega=Omega, inc=inc, T=T)

    # Add K2-19c
    m = row.masse3 * c.M_earth / c.M_sun # Solar masses
    per = row.per3
    a = (per**2 * sim.G * (Mstar + m) / 4 / np.pi**2)**(1/3.)
    e = row.e3
    omega = row.w3
    Omega = row.Omegarad3
    inc = row.incrad3
    T = row.tp3 # Time of periastron
    sim.add(m=m, a=a, e=e, omega=omega, Omega=Omega, inc=inc, T=T)
    
    sim.move_to_com()

    return sim

def plot_simulation(irow):
    sim = setupSimulation(irow)
    sim.status()
    sim.exit_min_distance = 0.01
    sim.automateSimulationArchive("archive.bin",interval=10,deletefile=True)
    tstart = sim.t
    tstop = 20000
    Noutputs = 1000

    times = np.linspace(tstart, tstop, Noutputs)
    a1 = np.zeros(Noutputs)
    e1 = np.zeros(Noutputs)
    M1 = np.zeros(Noutputs)
    w1 = np.zeros(Noutputs)
    a2 = np.zeros(Noutputs)
    e2 = np.zeros(Noutputs)
    M2 = np.zeros(Noutputs)
    w2 = np.zeros(Noutputs)

    distances = np.zeros(Noutputs)
    ps = sim.particles # ps is now an array of pointers. It will update as the simulation runs.

    try:
        for i,time in enumerate(times):
            sim.integrate(time, exact_finish_time=0)
            a1[i] = ps[1].a
            e1[i] = ps[1].e
            M1[i] = ps[1].M
            w1[i] = ps[1].omega
            a2[i] = ps[2].a
            e2[i] = ps[2].e
            M2[i] = ps[2].M
            w2[i] = ps[2].omega

            dp = ps[1] - ps[2]   # Calculates the coponentwise difference between particles
            distances[i] = np.sqrt(dp.x*dp.x+dp.y*dp.y+dp.z*dp.z)
    except rebound.Encounter as error:
        print(error)


    dw = w2 - w1
    dw = np.arctan2(np.sin(dw), np.cos(dw))
    fig1 = plt.figure(figsize=(15,5))

    ax = plt.subplot(131)
    ax.set_xlabel("time [days]")
    ax.set_ylabel("e1 (black), e2 (orange)")
    plt.plot(times, e1);
    plt.plot(times, e2);
    plt.ylim(0,0.3)

    ax = plt.subplot(132)
    ax.set_xlabel("$\Delta \omega$ [deg]")
    ax.set_ylabel("e1 (black), e2 (orange)")
    plt.plot(np.rad2deg(dw), e1,'.',ms=2);
    plt.plot(np.rad2deg(dw), e2,'.',ms=2);
    plt.xlim(-180,180)
    plt.ylim(0,0.3)

    ax = plt.subplot(133)
    ax.set_xlabel("time [days]")
    ax.set_ylabel("M1 - M2 [AU]")
    x = 3 * M2
    y = 2 * M1
    Lambda = np.arctan2(np.sin(x-y), np.cos(x-y))
    plt.plot(times,Lambda)
    plt.ylabel('3 * M2 - 2 * M1')

    fig2, axL = plt.subplots(ncols=2,figsize=(9,4),sharey=True)
    ax = axL[0]
    sa = rebound.SimulationArchive("archive.bin")
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



