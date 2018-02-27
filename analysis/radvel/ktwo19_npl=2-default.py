import radvel
import ktwo19.io
import numpy as np

starname = 'epic201505350'
nplanets = 2    # number of planets in the system
instnames = ['j']    # list of instrument names. Can be whatever you like but should match 'tel' column in the input file.
ntels = len(instnames)       # number of instruments with unique velocity zero-points
fitting_basis = 'per tc secosw sesinw k'    # Fitting basis, see radvel.basis.BASIS_NAMES for available basis names
bjd0 = 2454833
planet_letters = {1: 'b', 2:'c', 3:'d'}

aparams = radvel.Parameters(nplanets, basis='per tc secosw sesinw k') 

aparams['per1'] = radvel.Parameter(value=7.919520)
aparams['tc1'] = radvel.Parameter(value=1980.38319) # time of inferior conjunction of 1st planet
aparams['secosw1'] = radvel.Parameter(value=0.0)    
aparams['sesinw1'] = radvel.Parameter(value=0.0)    
aparams['k1'] = radvel.Parameter(value=10.0)         # velocity semi-amplitude for 1st planet

aparams['per2'] = radvel.Parameter(value=11.907244)
aparams['tc2'] = radvel.Parameter(value=1984.27545)    # time of inferior conjunction of 1st planet
aparams['secosw2'] = radvel.Parameter(value=0.0)    
aparams['sesinw2'] = radvel.Parameter(value=0.0)    
aparams['k2'] = radvel.Parameter(value=10.0)         # velocity semi-amplitude for 1st planet

aparams['dvdt'] = radvel.Parameter(value=0.0)        # slope
aparams['curv'] = radvel.Parameter(value=0.0)         # curvature
aparams['gamma_j'] = radvel.Parameter(1.0)      # velocity zero-point for hires_rj
aparams['jit_j'] = radvel.Parameter(value=2.6)        # jitter for hires_rj

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

