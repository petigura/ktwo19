import os
from cStringIO import StringIO

import numpy as np
import pandas as pd

import re
import ttv.lithwick
from cStringIO import StringIO as sio
import ktwo19.photometry
DATADIR = os.path.join(os.path.dirname(__file__),'../data/')
def load_table(table, cache=0, cachefn='load_table_cache.hdf', verbose=False):
    """Load tables

    Args:
        table (str): name of table. must be one of


        cache (Optional[int]): whether or not to use the cache
            - 0: don't use the cache recreate all files
            - 1: read from cache
            - 2: write tables to cache

    Returns:
        pandas.DataFrame: table

    """
    if cache==1:
        try:
            df = pd.read_hdf(cachefn,table)
            print "read table {} from {}".format(table,cachefn)
            return df
        except IOError:
            print "Could not find cache file: %s" % cachefn
            print "Building cache..."
            cache=2
        except KeyError:
            print "Cache not built for table: %s" % table
            print "Building cache..."
            cache=2

    if cache==2:
        df = load_table(table, cache=False)
        print "writing table {} to cache".format(table)
        df.to_hdf(cachefn,table)
        return df

    elif table=='ephem-sinukoff16':
        fn = 'data.xlsx'
        fn = os.path.join(DATADIR, fn)
        df = pd.read_excel(fn,sheet_name=table,usecols=[0,1],index_col=0,squeeze=True)
        for i in range(1,4):
            df['tc%i' % i] += (2456000 - 2454833) # Kepler epoch

    elif table=='ktwo-everest':
        df = ktwo19.photometry._everest()

    return df

def load_djh():
    s = """\
    tc               tc_err
    2457898.54880    0.00459 
    2457906.46919    0.00569 
    2457914.39057    0.00674 
    2457922.31341    0.00776 
    2457930.23357    0.00881 
    2457938.15467    0.00981 
    2457946.07734    0.01071 
    2457953.99728    0.01167 
    2457961.91812    0.01261 
    """
    s = pd.read_table(sio(s),sep='\s*')
    s['tc'] -=bjd0
    s['i_epoch'] = s.index
    s['i_epoch'] += 137
    s['i_planet'] = 1
    times_djh = s.copy()


    s = """\
    tc               tc_err
    2457900.04646    0.02629 
    2457911.94202    0.03194 
    2457923.83491    0.03837 
    2457935.73184    0.04370 
    2457947.62592    0.04959 
    2457959.52420    0.05443 
    2457971.41940    0.05959 
    """
    s = pd.read_table(sio(s),sep='\s*')
    s['tc'] -=bjd0
    s['i_epoch'] = s.index
    s['i_epoch'] += 91
    s['i_planet'] = 2

    times_djh = times_djh.append(s)
    times_djh.index = times_djh.i_planet
    return times_djh
