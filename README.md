# Estimation of multiple knockoff procedures
This package provides implementations of several multiple knockoff aggegration schemes in `R`. 

## Description
Multiple knockoff procedures run the knockoff filter multiple times, each time with a different knockoff matrix, and then aggregate the results in a way such that FDR control is (hopefully) still retained while reducing the variability from the probabilistic knockoff construction. The multiple knockoff filters also aim to increase the statistical power of the aggregated selection set.
The package implements the following three aggregation procedures for multiple knockoffs:
- Union knockoffs by Xie and Lederer (2021). This method runs multiple knockoff filters with different nominal levels whose sum equals the nominal level of our desired FDR control.
- p-value knockoffs by Nguyen et al. (2020). This aggregation scheme derives for each variable a p-value which is calculated from the multiple score statistic of choice. Then, the method applies either Benjamini-Hochberg or Benjamini-Yiekutieli to obtain the final selection set with FDR control.
- ADAGES by Gui (2020). This procedure aggregates multiple selection procedures that have FDR control at q respectively by finding an adaptive threshold integer. Then, ADAGES selects all variables that occur at least as often as the magnitude of the threshold across all selection sets.

## Installation
The package `multiknockoffs` can be directly installed in R with the devtools package by typing the following commands:
```
devtools::install_github("cKarypidis/multiknockoffs")
```

## Examples
Since this package is part of a Master's thesis, the user can find examples within that work or in the regarding help menu of the package in the `R` console. Similar to the `knockoff` package, which implements the fixed-X and model-X knockoff filter, the user can either run the whole procedure by one function or each step manually. In the following we will list the regarding sections of the Master's thesis that contain examples.
- Section 5.1.2 to run union knockoffs as whole procedure.
- Section 5.2.2 to run p-value knockoffs as whole procedure.
- Section 5.3.2 to run ADAGES as whole procedure.
- Section A.6.2 to run each step of the procedures manually. 

A vignette is in progress.

## Resources

Gui, Y. (2020). ADAGES. Proceedings of the 2020 ACM-IMS on Foundations of Data
Science Conference. ACM.
http://dx.doi.org/10.1145/3412815.3416881

Nguyen, T.-B., J.-A. Chevalier, B. Thirion, and S. Arlot (2020). Aggregation of Multiple
Knockoﬀs. Proceedings of the 37th International Conference on Machine Learning. Ed.
by H. D. III and A. Singh. Vol. 119. Proceedings of Machine Learning Research. PMLR,
7283–7293.
https://arxiv.org/abs/2002.09269?context=stat.ML

Xie, F. and J. Lederer (2021). Aggregating Knockoﬀs for False Discovery Rate Control
with an Application to Gut Microbiome Data. Entropy 23(2), 230. 
https://doi.org/10.3390/e23020230
