# Estimation of multiple knockoff procedures
This package provides implementations of several multiple knockoff aggegration schemes in R. 

## Description
Multiple knockoff procedures run the knockoff filter multiple times, each time with a different knockoff matrix, and then aggregate the results in a way such that FDR control is (hopefully) still retained while reducing the variability from the probabilistic knockoff construction. The multiple knockoff filters also aim to increase the statistical power of the aggregated selection set.
The package implements the following three aggregation procedures for multiple knockoffs:
- Union knockoffs by Xie and Lederer (2021). This method runs multiple knockoff filters with different nominal levels whose sum equals the nominal level at wich we want FDR control.
- p-value knockoffs by Nguyen et al. (2020). This aggregation scheme derives for each variable a p-value which is calculated from the multiple score statistic of choice. Then, the method applies either Benjamini-Hochberg or Benjamini-Yiekutieli to obtain a final selection set with FDR control.
- ADAGES by Gui (2020). This procedure runs aggregates multiple selection procedures that have FDR control at q respectively by finding an adaptive trheshold integer. Then, ADAGES selects all variables that occur at least as often as the magnitude of the threshold across all selection sets.

## Installation
The package `multiknockoffs` can be directly installed in R with the devtools package by typing the following commands:
```
devtools::install_github("cKarypidis/multiknockoffs")
```

## Examples


## Resources

