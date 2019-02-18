# K2-19 TTV+RV project

## Dependencies

The analysis doesn't work with the newest version of the code.

```
conda create -n ktwo19 python=2.7
source activate ktwo19
pip install radvel==1.2
pip install jinja2==2.0
conda install seaborn
conda install xarray
pip install celerite
pip install ChainConsumer
conda install matplotlib=2.0.2 
pip install xlrd
```

Run in the ktwo19 conda environment on Erik's laptop

## Receipe 

### 1. Gather data

Hand Entered data lives in `data/data.xlsx` where the sheets correspond to

- transit-times: list of transit times 
- ephem-sinukoff16: Ephemeris given in Sinukoff et al. (2016).  
- stellar-sme: stellar parameters from Brewer et al. (2016)

Get RV data with 

```
cp ~/.cpsutils/mir3/vel/vstepic201505350.dat data/
```

### 2. Derive stellar parameters

Update with Gaia values by following instructions in ``0_Stellar-Parameters.ipynb``

### 3. Train the GP using the photometry

```
1_Train-GP-with-photometry
```

### 4. Perform Keplerian RV fittng MCMC use

```
2_Model-Comparison-BIC
```

which will perform some radvel fits with models of varying complexity. Run radvel fits and MCMC. This script runs models of various complexities. Look at the reports.

NOTE: currently, this doesn't work because of radvel version
issues. If I need to run this again, just migrate the code to the
current version of radvel.

### 4. Prepare data for photodynamical modeling

```
run_ktwo19.py photdyn-prepare 
```

This generates the K2 photometry, which is stored

- `analysis/photodyn/photometry-ktwo.tsv` 

and the RVs which are stored

- `analysis/photodyn/rv.tsv`
- `analysis/photodyn/rv-trend-removed.tsv`

`rv.tsv` is the raw RV dataset. For `rv-trend-removed.tsv` I've subtracted off the MAP gamma a dvdt value.

- time:  BJD - 2454833
- mnvel: the mean RV
- errvel: measurement uncertainty

NOTE: the analysis did not model the RVs directly.

### 5 Run photodynamical code

Follow the instructions in `XX_Figure-Photodyn-Bestfit.ipynb`

### 5 Create the photodymm-bestfit figure

Follow the instructions in `XX_Figure-Photodyn-Bestfit.ipynb`

### 6 Create the TTV Fit figure

Follow the instructions here:

XX_Figure_TTVFit.ipynb
