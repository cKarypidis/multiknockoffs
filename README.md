# Estimation of multiple knockoff procedures
This package provides implementations of several multiple knockoff aggegration schemes in R. 

## Description
Multiple knockoff procedures run the knockoff filter multiple times, each time with a different knockoff matrix, and then aggregate the results in a way such that FDR control is (hopefully) still retained while reducing the variability from the probabilistic knockoff construction. The multiple knockoff filters also aim to increase the statistical power of the aggregated selection set.
The package implements the following three aggregation procedures for multiple knockoffs:
- Union-knockoffs by ... This method runs multiple knockoff filters with squence and ... 
- p-value knockoffs by 
- ADAGES by. 

## Installation
The package `multiknockoffs` can be directly installed in R with the devtools package by typing the following commands:
```
devtools::install_github("cKarypidis/multiknockoffs")
```

## Examples


## Resources

