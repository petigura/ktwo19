# K2-19 TTV+RV project

## Data

Hand entered data lives in

```
data/data.xlsx
```

where the sheets correspond to

- transit-times: list of transit times 
- ephem-sinukoff16: Ephemeris given in Sinukoff et al. (2016).  
- stellar-sme: stellar parameters from Brewer et al. (2016)

## Derive stellar parameters

Update with Gaia values

0_Stellar-Parameters.ipynb

cp ~/.cpsutils/mir3/vel/vstepic201505350.dat data/


### Perform Keplerian RV fitting

Run radvel fits and MCMC. This script runs models of various complexities. Look at the reports

```
run_ktwo24.py fit-rv-keplerian
cp ktwo24_npl\=3-ccc/ktwo24_npl\=3-ccc_rv_multipanel.pdf ../../paper/fig_rv-keplerian-circ.pdf 
```

Record the BIC for the in `paper/val_hand.tex`






## Notes for Sean

### Photometry

K2 Photometry (Everest)

```
analysis/photodyn/photometry-ktwo.tsv
```

### RVs

I've included two RV timeseries. 

```
analysis/photodyn/rv.tsv
analysis/photodyn/rv-trend-removed.tsv	
```

`rv.tsv` is the raw RV dataset. For `rv-trend-removed.tsv` I've subtracted off the MAP gamma a dvdt value.

- time:  BJD - 2454833
- mnvel: the mean RV
- errvel: measurement uncertainty


## Old notes

## Run MCMC  

```
$ bin/mcmc_mnest.py <method> <id> --livepoints=300
```

Currently, I'm getting good fits with 

method == 'ttvfast-npar9'

