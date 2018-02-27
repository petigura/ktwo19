# Derive stellar parameters

0_Stellar-Parameters.ipynb

cp ~/.cpsutils/mir3/vel/vstepic201505350.dat data/

# Run MCMC  

```
$ bin/mcmc_mnest.py <method> <id> --livepoints=300
```

Currently, I'm getting good fits with 

method == 'ttvfast-npar9'

