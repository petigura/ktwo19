import ktwo19.io
import everest 
import numpy as np
from matplotlib.pylab import *
from cpsutils.pdplus import LittleEndian as LE
import pandas as pd

def _everest():
    """
    Load and prepare everest photometry
    """
    ephem = ktwo19.io.load_table('ephem-sinukoff16')
    star = everest.Everest(201505350)

    # Mask out transits as recommended by Luger
    pad = 1.0
    star.mask_planet(ephem.tc1, ephem.per1, dur = (ephem.dur1 +pad) / 24 )
    star.mask_planet(ephem.tc2, ephem.per2, dur = (ephem.dur2 +pad) / 24 )
    star.mask_planet(ephem.tc3, ephem.per3, dur = (ephem.dur3 +pad) / 24 )
    star.compute()

    # Clip out outliers
    mask = np.zeros_like(star.time)
    mask[star.outmask] = 1
    mask[star.badmask] = 1
    mask[star.transitmask] = 0
    fm = ma.masked_array(star.flux,mask )
    fm.data[:] = fm.data / ma.median(fm)
    fm.mask = fm.mask | (fm < 0.98)
    plot(star.time,fm,'.')

    df = dict(t=star.time,f=fm.data,mask=fm.mask)
    df = pd.DataFrame(df)
    df = LE(df.to_records(index=False))
    df = df['t f mask'.split()]
    df = pd.DataFrame(df)
    return df
