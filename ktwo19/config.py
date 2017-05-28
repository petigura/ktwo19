DATADIR='/Users/petigura/Research/K2-24_TTV+RV/data/transit-times/'
bjd0=2454833

from astropy import constants as c 
from astropy import units as u
from numpy import exp

M_earth = (c.M_earth / c.M_sun).value 
AU_per_day = (1.0 * u.AU/u.d ).to(u.m/u.s).value

