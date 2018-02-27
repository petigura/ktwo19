"""
Three planets with circular orbits
"""
import os
import numpy as np
import radvel

import imp
from ktwo19.config import DATADIR
default = os.path.join(DATADIR,'../analysis/radvel/','ktwo19_npl=2-default.py')
P = imp.load_source('test', default)
for d in dir(P):
    if d[0]!='_':
        exec("{}  = P.{}".format(d,d))

params['per1'].vary = False # period of 1st planet
params['tc1'].vary = False 
params['secosw1'].vary = False
params['sesinw1'].vary = False
params['k1'].vary = True
params['per2'].vary = False
params['tc2'].vary = False
params['secosw2'].vary = False
params['sesinw2'].vary = False
params['k2'].vary = True
params['dvdt'].vary = True
params['curv'].vary= False
params['gamma_j'].vary = True
params['jit_j'].vary = True
