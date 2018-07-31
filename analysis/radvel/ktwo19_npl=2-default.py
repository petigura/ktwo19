import radvel
import ktwo19.io
import numpy as np

starname = 'epic201505350'
nplanets = 2    # number of planets in the system
instnames = ['j']    # list of instrument names. Can be whatever you like but should match 'tel' column in the input file.
ntels = len(instnames)       # number of instruments with unique velocity zero-points
fitting_basis = 'per tc secosw sesinw k'    # Fitting basis, see radvel.basis.BASIS_NAMES for available basis names
bjd0 = 2454833
planet_letters = {1: 'b', 2:'c'}

# Setup default values
aparams = radvel.Parameters(nplanets,basis='per tc secosw sesinw k')
aparams['per1'] = radvel.Parameter(value=7.920751,vary=False)
aparams['tc1'] = radvel.Parameter(value=1980.379057,vary=False)
aparams['sesinw1'] = radvel.Parameter(value=0.,vary=False) # fix eccentricity = 0
aparams['secosw1'] = radvel.Parameter(value=0.,vary=False)
aparams['k1'] = radvel.Parameter(value=10.0)
aparams['per2'] = radvel.Parameter(value=11.898036,vary=False)
aparams['tc2'] = radvel.Parameter(value=1984.304530,vary=False)
aparams['sesinw2'] = radvel.Parameter(value=0.,vary=False) # fix eccentricity = 0
aparams['secosw2'] = radvel.Parameter(value=0.,vary=False)
aparams['k2'] = radvel.Parameter(value=10.0)
aparams['dvdt'] = radvel.Parameter(value=0.0,vary=False)
aparams['curv'] = radvel.Parameter(value=0.0,vary=False)

params = aparams.basis.to_any_basis(aparams,fitting_basis)
# Load radial velocity data, in this example the data is contained in an hdf file,
# the resulting dataframe or must have 'time', 'mnvel', 'errvel', and 'tel' keys
# the velocities are expected to be in m/s
data = ktwo19.io.load_table('rv')
# abscissa for slope and curvature terms (should be near mid-point of time baseline)
time_base = int(data.time.mean())

# optional argument that can contain planet radii,
# used for computing densities. Values should be given
# in units of Earth radii

star = ktwo19.io.load_table('stellar')
stellar = dict(mstar=star.smass, mstar_err= star.smass_err)

# Note that these are temp values for planet radius in Earth-radii
planet = {
    'rp1':6.98,
    'rp_err1':0.33,
    'rp2':4.16,
    'rp_err2':0.33,
}


# Define prior shapes and widths here.
priors = [
    radvel.prior.EccentricityPrior( nplanets ), # Keeps eccentricity < 1
    radvel.prior.HardBounds('jit_j', 0.0, 10.0)
]

