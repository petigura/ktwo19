"""Hey Sean.

I did my best to convert the conversion code over into python. There
are a couple functions I'm missing.

- getlambda
- rotatekep
- getf

Would you also mind having a look at my pyhton code. I left some
questions and clafications as comments.

"""

def dsetup2(p, npl, epoch):
    """

    Args:
        p (array): vector that has for each planet. Array has length 8*np + 5 
            where:

            p[8*ipl + 0] : period (days)
            p[8*ipl + 1] : T0 (days)
            p[8*ipl + 2] : sqrt(e) cos (omega) 
            p[8*ipl + 3] : sqrt(e) sin (omega) 
            p[8*ipl + 4] : inclination (deg)
            p[8*ipl + 5] : Omega (deg)
            p[8*ipl + 6] : mass (jupiter masses)
            p[8*ipl + 7] : Rp/Rstar
            
            p[npl*pperplan+0] - stellar mass (solar masses)
            p[npl*pperplan+1] - stellar radius (solar radii)
            p[npl*pperplan+2] - c1 limb-darkening
            p[npl*pperplan+3] - c2
            p[npl*pperplan+4] - dillution

        epoch (float): osculating elements defined at this epoch

    Returns:

        array: (x, y, z, vx, vy, vz) values in a stellar-centric coordinate 
            system of each planet. 

    """

    # const int sofd = SOFD # Size of double? (don't need in Python)
    # const int sofds = SOFDS; # Size of double signed?
  
    pperplan = 8 # parameters per planet
    pstar = 5 # parameters per star

    ms = p[npl*pperplan + 0]
    rstar = p[npl*pperplan +1 ] # not needed in this function
    c1 = p[npl*pperplan + 2] # not needed in this function
    c2 = p[npl*pperplan + 3] # not needed in this function
    dilute = p[npl * pperplan+4] # not needed in this function
    
    bigg = 1.0e0 # Newton's constant # Not used in this function
    ghere = 2.9591220363e-4 # G # 2.9591220363e-4; Msun-AU-day units
    jos =  9.545e-4 # M_jup / M_sol
    
    mp = np.zeros(npl) # array for mass of planets 
    mpjup = np.zeros(npl) # array for mass of planets
    msys = np.zeros(npl+1) # mass of system
    msys[0] = ms # mass of star is first element in msys

    a = np.zeros(np) # semi-major axis [AU]
    e = np.zeros(npl) # eccentricity
    inc = np.zeros(npl) # inclination [rad]
    bo = np.zeros(npl) # longitude of ascending node [rad]
    lo = np.zeros(npl) # argument of peri [rad]
    lamb = np.zeros(npl) # mean longitude [rad]
    f = np.zeros(npl) # true anomaly [rad]
    M_PI = np.pi

    for i in range(npl):
        per = p[i*pperplan+0]
        T0 = p[i*pperplan+1]
        secosw = p[i*pperplan+2]
        sesinw = p[i*pperplan+3]

        e[i] = secosw**2 + sesinw**2

        inc[i] = p[i*pperplan+4] * M_PI/180.0 # converting to radians
        bo[i] = p[i*pperplan+5] * M_PI/180.0 # converting to radians
        lo[i] = np.arctan2(sesinw,secosw) # argument of peri
        mp[i]= p[i*pperplan+6] # mass of planet

        # Convert mp [Mjup] to [Msol] (mpjup is a confusing variable name)
        mpjup[i] = mp[i] * jos 
        msys[i+1] = msys[i] + mpjup[i] # mass of body, and all internal bodies

        # Kepler's third law
        a[i] = (ghere * msys[i+1] * per**2.0 / 4.0 / M_PI**2)**(1.0/3.0)
        pomega = bo[i] + lo[i] # longitude of pericenter [rad]

        # Is this the mean longitude at inferior conjunction? I don't
        # have the getlambda function. Transit will occur when 
        # f + lo = pi/2. Therefore, getlambda needs f, e, and pomega, and
        # Involves solving Kepler's equation.
        lamb0 = getlambda( (M_PI/2.0 - lo[i]), e[i], pomega)
        m0 = lamb0 - pomega; # Mean anomaly at time of transit
        me = m0 + 2 * M_PI * (epoch - T0)/per # Mean anomaly at epoch
        mepomega = me + pomega # Mean longitude at epoch
        lambepoint = pushpi(&mepomega,1) # Mean longitude at epoch?
        lambe = lamepoint[0] # not sure what this does
        f[i] = getf(lambe, e[i], pomega) # I need the function getf

    state = np.zeros(npl * pperplan)
    for i in range(npl):
        stateplan = keptostate(
            a[i],e[i],inc[i],lo[i],bo[i],f[i],ghere*msys[i+1]
        )
        state[i*pperplan:i*pperplan+pperplan] = stateplan

    state = -1.0 * state # Not sure why it needs to be inverted

def keptostate(a, e, i, lo, bo, f, m):
    """

    Args:
        a (float): semi-major axis
        e (float): eccentricity
        i (float): inclination [rad]
        lo (float): argument of peri "litle omega" [rad]
        bo (float): longitude of ascending node "big omega" [rad]
        f (float): true anomaly [rad]
        m (float): mass times G 

    Returns:
        array: x, y, z, vx, vy, vz 

    """
    mass = m
    r = a*(1.0 - e**2)/(1.0 + e*np.cos(f)) # radial distance
    x0 = r * np.cos(f)
    y0 = r * np.sin(f)
    vx0 = -np.sqrt(mass/( a * (1.0-e**2))) * np.sin(f)
    vy0 = np.sqrt(mass/(a * (1.0 - e**2))) * (e+np.cos(f))

    statearr = np.zeros(6)
    state1arr = rotatekep(x0, y0, i, lo, bo)
    state2arr = rotatekep(vx0, vy0, i, lo, bo)
    statearr = np.hstack([state1arr, state2arr])
    return statearr



