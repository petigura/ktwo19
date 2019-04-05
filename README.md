# K2-19 TTV+RV project

## Dependencies

The analysis doesn't work with the newest version of radvel. I had to create a conda environment `ktwo19-2` and install the following packages

```
conda create -n ktwo19-2 python=2.7
source activate ktwo19
pip install radvel==1.2
pip install jinja2==2.0
conda install seaborn
conda install xarray
pip install celerite
pip install ChainConsumer
conda install matplotlib=2.0.2 
pip install xlrd
conda install pytables
```

also needed to install rebound used ef325d3dd69445c25aeeb4a6b020ac18cde8ac4b

There was also an issue with the macosx backend, so I installed a
matplotlibrc file and had it use tkagg


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

### 4. Perform Keplerian RV fittng

First compare models using 

```
jupyter notebook 2_Model-Comparison-BIC.ipynb
```

Model 5. with the circular planet b and the GP noise model is
preferred. However, run Model 6 with the eccentric planet b to compute
an upper bound on eccentricity.

Here's what I found:

```
                  mean   std
k1              11.369 1.829
k2               0.774 1.717
k3               0.194 1.237
secosw1         -0.133 0.163
sesinw1          0.144 0.224
dvdt             0.016 0.007
gp_amp           7.100 2.396
gp_explength    35.149 4.645
gp_per          20.339 0.315
gp_perlength     0.470 0.065
gamma_j         -2.265 2.167
jit_j            4.535 1.376
lnprobability -184.327 2.623
```

Which is within about 1.5 sigma of the photodynamical model and gives
a 95\% upper limit of 0.27

The final RV fit, which goes in the paper is measured through the 

```
ktwo19.keplerian.mcmc() method
```

### 5. Prepare data for photodynamical modeling

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

Set up runs on cadence

```
/mir3/petigura/code/Phodymm/
```

Run the code (documentation needed) 

edit `analysis/photodyn/runs.txt` with the run name

```
bin/pull_photdym.sh
```

### 5 Create the phodymm-bestfit figure

Follow the instructions in `XX_Figure-Photodyn-Bestfit.ipynb`

### 6 Create the TTV Fit figure

Follow the instructions here:

XX_Figure_TTVFit.ipynb

### 7 Run the N-body simulations

run_ktwo19.py  create-nbody-batch | grep run > nbody.tot
parallel :::: nbody.tot

