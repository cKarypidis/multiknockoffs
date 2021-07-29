# Estimation of multiple knockoff procedures
This package provides implementations of several multiple knockoff aggegration schemes in R. 

## Description
The knockoff filter (Zitate) is a procedure that achieves finite sample FDR control for almost all model classes. The probabilistic construction of (model-X) knockoffs introduces randomness. Running the knockoff algorithm multiple times on the same data results in different knockoff matrices and probably overlapping but slightly different selection sets. These differences will yield to fluctuating power and FDP values across different runs which make inference more irreproducible. Hence, a knockoff filter with a more stable selection of variable is desirable. Multiple knockoff procedures consist of running the knockoff filter multiple times, each time with a different knockoff matrix, and then aggregate the results in a way such that FDR control is still retained while reducing the variability from the knockoff construction. The multiple knockoff aggregation filters do not lead to a stable selection but they can also increase the statistical power.

## Installation
The package `multiknockoffs` can be directly installed in R with the devtools package by typing the following commands:
```
devtools::install_github("cKarypidis/multiknockoffs")
```

## Examples


## Resources

