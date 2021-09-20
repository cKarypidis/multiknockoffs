# Estimation of multiple knockoff procedures
The package `multiknockoffs` provides implementations of several multiple knockoff aggegration schemes in `R`. 

## Description
The knockoff filter (Barber and Candès (2015); Candès et al. (2018)) is a modern and powerful algorithm to control the false discovery rate (FDR) for a variety of different model classes, including complex machine learning models such as neural networks, boosting, random forests or high-dimensional linear penalization methods. The procedure constructs fake features (a knockoff matrix), that mimic certain correlation properties of the original variables. Since knockoffs behave similar to the original
features but are known to be artificial null variables, they serve as a negative control
group.

Multiple knockoff procedures run the knockoff filter multiple times, each time with a different generated knockoff matrix, and then aggregate the results to
to improve the power while (approximately) retaining empirical FDR control of the aggregated selection set. The multiple knockoff filters also aim to reduce the variability of  power and empirical FDR values resulting from the probabilistic knockoff construction. The package implements the following three aggregation procedures for multiple knockoffs:
- Union knockoffs by Xie and Lederer (2021). This method runs multiple knockoff filters with different nominal levels whose sum equals the nominal level of our desired FDR control.
- p-value knockoffs by Nguyen et al. (2020). This aggregation scheme derives for each variable a p-value which is calculated from the multiple score statistic of choice. Then, the method applies either Benjamini-Hochberg or Benjamini-Yiekutieli to obtain the final selection set with FDR control.
- ADAGES by Gui (2020). This procedure aggregates multiple selection procedures that have FDR control at q respectively by finding an adaptive threshold integer. Then, ADAGES selects all variables that occur at least as often as the magnitude of the threshold across all selection sets.

Since the package has resulted complementary to my [Master's thesis](documentations/Masterthesis.pdf), the work contains compact explanations of each method. The user can have a closer look at Section 5.1 for union knockoffs, Section 5.2 for p-value knockoffs, Section 5.3 for ADAGES and Section 5.4 for a simulation study that compares the three methods regarding their power and empirical FDR values. Finally, Appendix B of the work presents some additional implementation details.

## Installation
The package `multiknockoffs` can be directly installed in `R` with the devtools package by typing the following commands:
```
devtools::install_github("cKarypidis/multiknockoffs")
```

## Examples
Similar to the `knockoff` package, which implements the fixed-X and model-X knockoff filter, the user can either run the whole procedure by one function or each step manually.
In the following, we will present the application of each multiple knockoff method by the "all-in-one" approach. 

Assume we generate data according to
```
n <- 400; p <- 200; s_0 <- 30   #s_0 is the number of true signal variables where p is the total number of variables
amplitude <- 1; mu <- rep(0,p); rho <- 0.25
Sigma <- toeplitz(rho^(0:(p-1)))

X <- MASS::mvrnorm(n, mu, Sigma)
nonzero <- sample(p, s_0)
beta <- amplitude * (1:p %in% nonzero)
y <- X %*% beta + rnorm(n)
```
The three multiple knockoff procedures with default options can be applied respectively by

```
#Union knockoffs
res.uKO <- run.uKO(X, y, sets = TRUE)

#P-value knockoffs
res.pKO <- run.pKO(X, y, pvals = TRUE)

#ADAGES
res.ADAGES <- run.ADAGES(X, y, sets = TRUE)
```

For more in-depth explanations of the functions' arguments and outputs as well as the manual approach of each step, we refer to the **vignette** pdf document. Moreover, the user can find explanations and examples in the regarding help menu of the package in the `R` console. 

A vignette is in progress.

## Resources

**Original knockoffs:**

Barber, R. F. and E. J. Candès (2015). Controlling the false discovery rate via knockoffs.
The Annals of Statistics 43(5), 2055-2085.
http://dx.doi.org/10.1214/15-AOS1337

Candès, E., Y. Fan, L. Janson, and J. Lv (2018). Panning for gold: 'model-X' knockoffs for
high dimensional controlled variable selection. Journal of the Royal Statistical Society:
Series B (Statistical Methodology) 80(3), 551-577.
https://doi.org/10.1111/rssb.12265


**Multiple knockoffs:**

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
