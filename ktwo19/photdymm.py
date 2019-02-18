from __future__ import print_function
# DEMCMC analysis
# Usage:
# $ python demcmc_quick_analyze.py kepler36_longcadence.in [burnin]
# 
# burnin is an optional long integer of the number of generations to discard at the beginning of the demcmc output file
#

##
# additional uncertainty in stellar mass to include 
msigma = 0.0
##

import sys
import os
import pandas as pd
import numpy as np
import corner
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def read_phodymm(fname, demcmcfile, burnin):
    burnin //= 100
    f=open(fname)
    lines=f.readlines()

    def get_in_value(linenumber):
      line = lines[linenumber]
      line = line.split()
      return line[2]

    runname = get_in_value(8)
    #runname = runname[:-1]
#    demcmcfile='demcmc_'+runname+'.out'
    nchain = int(get_in_value(14))
    npl = int(get_in_value(10)) - 1
    parlist = lines[91+npl]
    npar = len(parlist.split())
    ndead=1



    # lines per output in demcmc_ file
    nper=npl+npar+ndead
    # string to save files with
    #savestr = demcmcfile.partition('demcmc_')[2]
    #savestr = savestr.partition('.out')[0]
    savestr = runname

    print(nper, npl, npar, ndead)

    # Read in file
    try:
      df1 = pd.read_csv(demcmcfile, sep='\t', header=None, skiprows=lambda x: x % nper >= npl, comment=';')
    except IOError:
      print("Error: This script must be run in the same directory as the 'demcmc' output file:")
      print("    "+demcmcfile)
    except Exception:
      print("The .in file was apparently parsed incorrectly, or your demcmc_NAME.out file has been corrupted.")
      print("   The .in file parsing has obtained the following values:")
      print("   N planets = %i" % npl)
      print("   N non-planet parameters (5 stellar + jitter + GPs) = %i" % npar)
      print("   N total rows per parameter set in the demcmc_NAME.out file = %i" % nper)
      print("")


    df2 = pd.read_csv(demcmcfile, sep='\t', header=None, skiprows=lambda x: (x % nper < npl or x % nper >= npl+npar), comment=';')


    # Number of parameters for each planet
    pperplan=df1.iloc[0].size

    # Restructure data into dataframe that has row indexes corresponding to generation
    #   and columns corresponding to different parameters
    allchs=[]
    for k in range(nchain):
      gens = [df1.iloc[k*npl+i::npl*nchain] for i in range(npl)]
      [i.reset_index(drop=True, inplace=True) for i in gens]
      chk = pd.concat(gens, axis=1, ignore_index=True) 

      gens_extra = [df2.iloc[k*npar+i::npar*nchain] for i in range(npar)]
      [i.reset_index(drop=True, inplace=True) for i in gens_extra]
      chk=pd.concat([chk]+gens_extra, axis=1, ignore_index=True)

      chk=pd.concat([chk]+[pd.DataFrame(k, index=np.arange(chk.shape[0]), columns=['Chain#'])], axis=1, ignore_index=True)

      allchs = allchs + [chk]

    allparam = pd.concat(allchs)


    ## Construct list of column names
    ##pparamlist=["Planet", "Period", "T0", "sqrt(e)*cosw", "sqrt(e)*sinw", "i", "Omega", "Mjup", "Rp/Rs"]
    ##if pperplan > 9:
    ##  pparamlist += ["bright", "c_1", "c_2"]
    ##sparamlist=["Ms", "Rs", "c_1", "c_2", "dilute"]
    ##extralist=["sigma"]
    ## Construct list of column names with fancy text
    pparamlist=["Planet", "Period", "T$_0,$", r"$\sqrt{e}\cos \omega$", r"$\sqrt{e}\sin \omega$", "i", r"$\Omega$", "M$_{jup,}$", "R$_p$/R$_s$"]
    if pperplan > 9:
      pparamlist += ["bright", "c$_1$", "c$_2$"]
    sparamlist=["M$_s$", "R$_s$", "c$_1$", "c$_2$", "dilute"]
    extralist=["$\sigma$"]

    fullplist=[]
    alphabet = 'bcdefghijklmnopqrstuvwxyz'
    for i in range(npl):
      fullplist += [pi+'$_'+alphabet[i]+'$' for pi in pparamlist]
    fullplist += sparamlist
    for i in range(npar - len(sparamlist)):
      fullplist += [a+"$_"+alphabet[i]+'$' for a in extralist] 
    fullplist += ["Chain#"]


    #allparam.to_csv('analysis_dir/out_4_all_preb.txt', sep='\t')
    # Remove burnin
    allparam.drop(range(burnin), inplace=True)

    allparam.columns = fullplist





    return allparam

def read_photdymm_format(_dir, tag):
    fname = '{}/{}.in'.format(_dir,tag)
    demcmcfname = '{}/demcmc_{}.out'.format(_dir,tag)
    df = ktwo19.photdymm.read_phodymm(fname, demcmcfname, 500)
    namemap = {}
    oldcols = ['Period','T$_0,$','$\sqrt{e}\cos \omega$','$\sqrt{e}\sin \omega$','i','$\Omega$','M$_{jup,}$','R$_p$/R$_s$']
    newcols = ['per','tc','secosw','sesinw','inc','Omega','Mjup','ror']

    for let, i  in zip(['$_b$','$_c$','$_d$'],['1','2','3']):
        for o,n in zip(oldcols,newcols):
            namemap[o+let] = n + i

    df = df.rename(columns=namemap)
    for let, i  in zip(['$_b$','$_c$','$_d$'],['1','2','3']):
        df['e'+i] = df['secosw'+i]**2 + df['sesinw'+i]**2
        df['ecosw'+i] = sqrt(df['e'+i]) * df['secosw'+i]
        df['esinw'+i] = sqrt(df['e'+i]) * df['sesinw'+i]
        df['masse' + i] = df['Mjup' + i] * 316 

