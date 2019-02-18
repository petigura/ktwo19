"""Hey Sean.

I did my best to convert the conversion code over into python. There
are a couple functions I'm missing.

- getlambda
- rotatekep
- getf

Would you also mind having a look at my pyhton code. I left some
questions and clafications as comments.

"""
import radvel._kepler
from  matplotlib.pylab import * 
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
  
    pperplan = 8 # parameters per planet
    pstar = 5 # parameters per star

    ms = p[npl*pperplan + 0]
    rstar = p[npl*pperplan +1 ] # not needed in this function
    c1 = p[npl*pperplan + 2] # not needed in this function
    c2 = p[npl*pperplan + 3] # not needed in this function
    dilute = p[npl * pperplan+4] # not needed in this function
    
    ghere = 2.9591220363e-4 # G # 2.9591220363e-4; Msun-AU-day units
    jos =  9.545e-4 # M_jup / M_sol
    
    mp = np.zeros(npl) # array for mass of planets 
    mpjup = np.zeros(npl) # array for mass of planets
    msys = np.zeros(npl+1) # mass of system
    msys[0] = ms # mass of star is first element in msys

    a = np.zeros(npl) # semi-major axis [AU]
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

        lamb0 = getlambda( (M_PI/2.0 - lo[i]), e[i], pomega) # lambda at the transit
        m0 = lamb0 - pomega; # Mean anomaly at time of transit
        me = m0 + 2 * M_PI * (epoch - T0)/per # Mean anomaly at epoch
        mepomega = me + pomega # Mean longitude at epoch
        lambe = pushpi(mepomega) # Mean longitude at epoch?
        f[i] = getf(lambe, e[i], pomega) 

    state = np.zeros(npl * 6)
    for i in range(npl):
        stateplan = keptostate(
            a[i],e[i],inc[i],lo[i],bo[i],f[i],ghere*msys[i+1]
        )
        state[i*6:i*6+6] = stateplan

    state = -1.0 * state # Inverted due to our choice of xyz orientation
    
    return state


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


def getlambda(f, e, pomega):
    bigE = 2.0 * np.arctan( np.sqrt((1.-e) / (1.+e)) * np.tan(f / 2.) )
    lam = pomega + bigE - e * np.sin(bigE)
    return lam

def getf(lam, e, pomega):
    bigM = lam - pomega
    bigE = bigM
    for i in range(20):
        bigE = bigM + e * np.sin(bigE)
    f = 2.0 * np.arctan( np.sqrt((1.+e) / (1.-e)) * np.tan(bigE / 2.) )
    return f

def pushpi(angle):
    while angle > np.pi:
        angle -= (2.*np.pi)
    while angle < -np.pi:
        angle += (2.*np.pi)
    return angle

def rotatekep(x, y, i, omega, bigO):
    z = 0.
    x1 = np.cos(omega) * x - np.sin(omega) * y
    y1 = np.sin(omega) * x + np.cos(omega) * y
    z1 = z
    x2 = x1
    y2 = np.cos(i) * y1 - np.sin(i) * z1
    z2 = np.sin(i) * y1 + np.cos(i) * z1
    x3 = np.cos(bigO) * x2 - np.sin(bigO) * y2
    y3 = np.sin(bigO) * x2 + np.cos(bigO) * y2
    z3 = z2
    return x3, y3, z3


def compute_true_anomaly(t, per, a, e, tp):
    """
    Compute true anomaly at specified time

    Args:
        t (float): time to compute true anomaly
        per (float): orbital period
        a (float): semimajor axis 
        e (float): eccentricity 
        tp (float): time of periastron passage

    Returns
        float: true anomaly
    
    """
    n = 2 * pi / per
    M = n * (t - tp) # compute mean anomaly
    M = np.array([M]) # needs array
    E = radvel._kepler.kepler_array(M,e) # solve kepler's equation
    E = E[0]
    f = 2 * np.arctan(np.sqrt((1+e)/(1-e)) * tan(E/2) ) # true anomally
    return f


def rotate_cartesian(cart, inc, Omega, omega):
    """
    cart: tuple (x,y,z,vx,vy,vz)
    """
    
    P1,P2,P3 = rotation_matricies(inc, Omega, omega)
    x,y,z,vx,vy,vz = cart

    pos0 = np.matrix(
        [[x],
         [y],
         [z]]
    )
    vel0 = np.matrix(
        [[vx],
         [vy],
         [vz]]
    )
    pos1 = P3 * P2 * P1 * pos0
    vel1 = P3 * P2 * P1 * vel0
    
    x1 = pos1[0,0]
    y1 = pos1[1,0]
    z1 = pos1[2,0]
    vx1 = vel1[0,0]
    vy1 = vel1[1,0]
    vz1 = vel1[2,0]
    return x1, y1, z1, vx1, vy1,vz1

def rotation_matricies(inc, Omega, omega):
    """
    Rotates coordinate by omega, inc, and Omega  MD p50
    """
    P1 = np.matrix([
        [cos(omega), -sin(omega), 0],
        [sin(omega), cos(omega), 0],
        [0, 0, 1],
    ])

    P2 = np.matrix([
        [1, 0, 0],
        [0, cos(inc), -sin(inc)],
        [0, sin(inc), cos(inc)],
    ])

    P3 = np.matrix([
        [cos(Omega), -sin(Omega), 0],
        [sin(Omega), cos(Omega), 0],
        [0,0,1]
    ])
    return P1, P2, P3

def keplerian_to_cartesian(per, a, e, inc, Omega, omega, f):
    """

    Confirmed that this works with worked example in MD p51
    """
    n = 2 * pi / per
    r = a * (1 - e**2) / (1 + e * cos(f))
    x = r * cos(f)
    y = r * sin(f)

    C = n * a / np.sqrt(1 - e**2) # constant 
    vx = - C * sin(f)
    vy = C * (e + cos(f))

    cart = x,y,0,vx,vy,0
    cart1 = rotate_cartesian(cart, inc, Omega, omega)
    return cart1

def cartesian_to_keplerian(cart, msys):
    x,y,z,vx,vy,vz = cart

    hx = y*vz - z*vy  # 2.129
    hy = z*vx - x*vz
    hz = x*vy - y*vx
    h = np.sqrt(hx**2 + hy**2 + hz**2)

    R = sqrt(x**2 + y**2 + z**2) # eqn 2.126
    V = sqrt(vx**2 + vy**2 + vz**2) # eqn 2.126

    RVR = x*vx + y*vy + z*vz # 2.218
    VR = sign(RVR) * sqrt(V**2 - h**2 / R**2)

    G = 2.9591220363e-4 # G # 2.9591220363e-4; Msun-AU-day units

    a = (2 / R - V**2 / (G * msys))**-1  # 2.134
    e = np.sqrt(1 - h**2 / G / msys/ a) # 2.135
    inc = arccos(hz / h)

    hsininc = h * sin(inc)
    # step 4
    if hz > 0:
        sinOmega = hx / hsininc
        cosOmega = -hy / hsininc
    else:
        sinOmega = -hx / hsininc
        cosOmega = +hy / hsininc

    if hsininc==0:
        print "inc == 0, so node is undefined"
        Omega = 0.0
    else:
        Omega = arctan2(sinOmega, cosOmega)
    
    # Step 5
    sinomegaf = z / R / sin(inc) # sin(omega + f)
    cosomegaf = 1 / cos(Omega) * (x / R + sin(Omega)*sinomegaf * cos(inc) ) 
    omegaf = arctan2(sinomegaf,cosomegaf)

    # Step 6
    sinf = a * (1 - e**2) / h / e * VR
    cosf = 1 / e * (a *(1 - e**2)/R - 1 )

    f = arctan2(sinf,cosf)
    omega = omegaf - f 
    return a, e, inc, Omega, omega, f 

