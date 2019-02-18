# All circular (preferred model with preferred bic)
radvel fit -s ktwo19_npl=2-cc.py 
radvel mcmc -s ktwo19_npl=2-cc.py 
radvel derive -s ktwo19_npl=2-cc.py 
radvel bic -t nplanets -s ktwo19_npl=2-cc.py 
radvel plot -t rv corner derived -s ktwo19_npl=2-cc.py 
radvel report -s ktwo19_npl=2-cc.py
