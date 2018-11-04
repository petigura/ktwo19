import rebound
import numpy as np
from astropy import constants as c
import radvel.orbit
import ktwo19.io
from matplotlib import pylab as plt

def setupSimulation(i):
    df = ktwo19.io.load_table('photodyn-samples',cache=1)
    row = df.iloc[i]
    sim = rebound.Simulation()
    sim.t = 2020.0 
    sim.G = 2.9591220363e-4 # in Msun-AU-day units
    sim.integrator = "ias15" 
    Mstar = 0.88
    sim.add(m=0.88) # K2-19 
    
    # Add K2-19b
    m = row.masse2 * c.M_earth / c.M_sun # Solar masses
    P = row.per2
    e = row.e2
    omega = row.w2
    Omega = row.Omegarad2
    inc = row.incrad2
    T = row.tp2 # Time of periastron
    sim.add(m=m, P=P, e=e, omega=omega, Omega=Omega, inc=inc, T=T)

    # Add K2-19c
    m = row.masse3 * c.M_earth / c.M_sun # Solar masses
    P = row.per3
    e = row.e3
    omega = row.w3
    Omega = row.Omegarad3
    inc = row.incrad3
    T = row.tp3 # Time of periastron
    sim.add(m=m, P=P, e=e, omega=omega, Omega=Omega, inc=inc, T=T)
    sim.move_to_com()

    return sim

def run_simulation(irow):
    sim = setupSimulation(irow)
    sim.status()
    sim.exit_min_distance = 0.01
    sim.automateSimulationArchive("archive.bin",interval=10,deletefile=True)
    tstart = sim.t
    tstop = 50000
    Noutputs = 1000

    times = np.linspace(tstart, tstop, Noutputs)
    a1 = np.zeros(Noutputs)
    e1 = np.zeros(Noutputs)
    M1 = np.zeros(Noutputs)
    n1 = np.zeros(Noutputs)
    w1 = np.zeros(Noutputs)
    a2 = np.zeros(Noutputs)
    e2 = np.zeros(Noutputs)
    M2 = np.zeros(Noutputs)
    n2 = np.zeros(Noutputs)
    w2 = np.zeros(Noutputs)

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
            a2[i] = ps[2].a
            e2[i] = ps[2].e
            M2[i] = ps[2].M  # mean anomaly
            n2[i] = ps[2].n  # mean motion
            w2[i] = ps[2].omega # argument
            
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
    #Lambda = np.arctan2(np.sin(x-y), np.cos(x-y))
    plt.plot(times,2*n1-3*n2)
    plt.ylabel(r'$2 n_1 - 3 n_2$ (rad/day) (black), $\dot{\omega}_1$ (rad/day) (red), ')
    plt.tight_layout(True)

    omegadot = (w1[1:] - w1[:-1]) / (times[1:] - times[:-1])
    Delta = n1/n2 * 2/3 - 1
    plt.plot(times,Delta)

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



