---
editor_options: 
  markdown: 
    wrap: 72
---

# multicenterCompRisk

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/multicenterCompRisk)](https://CRAN.R-project.org/package=multicenterCompRisk)
[![R-CMD-check](https://github.com/DELMASben/multicenterCompRisk/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DELMASben/multicenterCompRisk/actions/workflows/R-CMD-check.yaml)
[![License:
GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<!-- badges: end -->

## 📚 Overview

The goal of **`multicenterCompRisk`** is to provide tools for
**competing risks survival analysis** in **multicenter studies**,
accounting for **unobserved heterogeneity** (frailty) at the center
level as proposed in the following article:

Katsahian S, Resche-Rigon M, Chevret S, Porcher R. Analysing multicentre
competing risks data with a mixed proportional hazards model for the
subdistribution. *Stat Med*. 2006 Dec 30;25(24):4267-78. doi:
[10.1002/sim.2684](https://onlinelibrary.wiley.com/doi/10.1002/sim.2684).
PMID: [16960919.](https://pubmed.ncbi.nlm.nih.gov/16960919/)

This package is particularly useful for statisticians and
epidemiologists analyzing time-to-event data where individuals are
nested within centers (e.g., hospitals or clinics), and where multiple
causes of failure may occur.

It includes functions to: - Simulate clustered competing risks data, -
Fit cause-specific Cox models with shared frailty, following the method
presented by Sandrine Katsahian in "Analysing multicentre competing
risks data with a mixed proportional hazards model for the
subdistribution." - Fit cause-specific Cox models, - Fit Cox models with
shared frailty, - Fit Cox models.

The REML method is used when shared frailty is modeled; otherwise, the
ML method is applied.

## Results

We conducted the same simulation-based tests as in "Analysing
multicentre competing risks data with a mixed proportional hazards model
for the subdistribution.", and we obtained similar results:

| N   | K   | θ   | θ̂     | MSE(θ̂)  | γ       | γ̂        | MSE(γ̂)    |
|-----|-----|-----|-------|---------|---------|----------|-----------|
| 200 | 5   | 0.0 | 0.024 | 0.00244 | (0)     | (-0.068) | (0.23858) |
| 200 | 20  | 0.0 | 0.055 | 0.01027 | (0)     | (0.06)   | (0.34056) |
| 500 | 5   | 0.0 | 0.012 | 0.00062 | (0)     | (-0.001) | (0.10598) |
| 500 | 20  | 0.0 | 0.020 | 0.00137 | (0)     | (0.004)  | (0.0921)  |
| 200 | 5   | 0.0 | 0.036 | 0.00690 | (0.375) | (0.345)  | (0.23996) |
| 200 | 20  | 0.0 | 0.054 | 0.01221 | (0.375) | (0.381)  | (0.29773) |
| 500 | 5   | 0.0 | 0.012 | 0.00068 | (0.375) | (0.354)  | (0.10479) |
| 500 | 20  | 0.0 | 0.020 | 0.00161 | (0.375) | (0.42)   | (0.10242) |
| 200 | 5   | 0.1 | 0.130 | 0.03940 | (0)     | (-0.069) | (0.25903) |
| 200 | 20  | 0.1 | 0.127 | 0.02099 | (0)     | (-0.025) | (0.31313) |
| 500 | 5   | 0.1 | 0.108 | 0.01366 | (0)     | (-0.013) | (0.0985)  |
| 500 | 20  | 0.1 | 0.106 | 0.00754 | (0)     | (-0.056) | (0.10335) |
| 200 | 5   | 0.1 | 0.127 | 0.02148 | (0.375) | (0.372)  | (0.26376) |
| 200 | 20  | 0.1 | 0.139 | 0.02608 | (0.375) | (0.365)  | (0.27877) |
| 500 | 5   | 0.1 | 0.111 | 0.01158 | (0.375) | (0.378)  | (0.10493) |
| 500 | 20  | 0.1 | 0.099 | 0.00713 | (0.375) | (0.377)  | (0.10856) |
| 200 | 5   | 0.6 | 0.655 | 0.43517 | (0)     | (0.035)  | (0.22487) |
| 200 | 20  | 0.6 | 0.541 | 0.11904 | (0)     | (-0.009) | (0.35059) |
| 500 | 5   | 0.6 | 0.611 | 0.29234 | (0)     | (0.001)  | (0.10746) |
| 500 | 20  | 0.6 | 0.598 | 0.06516 | (0)     | (0.004)  | (0.10813) |
| 200 | 5   | 0.6 | 0.649 | 0.41894 | (0.375) | (0.398)  | (0.2656)  |
| 200 | 20  | 0.6 | 0.570 | 0.11168 | (0.375) | (0.378)  | (0.25024) |
| 500 | 5   | 0.6 | 0.629 | 0.30846 | (0.375) | (0.381)  | (0.10446) |
| 500 | 20  | 0.6 | 0.569 | 0.07018 | (0.375) | (0.38)   | (0.09621) |

## 📦 Installation

You can install the development version of `multicenterCompRisk` from
GitHub using either [`pak`](https://pak.r-lib.org/) or
[`devtools`](https://github.com/r-lib/devtools):

``` r
# Using pak (recommended)
# install.packages("pak")
pak::pak("DELMASben/multicenterCompRisk")

# Or using devtools
# install.packages("devtools")
devtools::install_github("DELMASben/multicenterCompRisk")
```

## 📌 Usage

### 🧪 Data simulation

You can generate simple datasets following the competing risks model
with center effect frailty with the code below :

``` r
library(multicenterCompRisk)

## SIZE AND CENTER REPARTITION
n_per_cluster=20 #sample size of each center
n_cluster=5 #number of centers
n=n_per_cluster * n_cluster #sample size
G=rep(1:n_cluster, each = n_per_cluster) #center repartition

## COVARIATES AND COMPETING FAILURE TIMES
n_cov=5 #number of covariables
Z=matrix(rnorm(n*n_cov,0,1),ncol = n_cov) #standardized Gaussian covariables 
prop=0.6 #proportion of failure 1
beta=c(1,1.2,0,-0.5,-0.2) #effects of covariates on failure cause 1
theta=0.6 #Variance of center effect frailty

## CENSORING
cens=TRUE #censoring in the simaluted data
pcens=0.25 #proportion of censoring
tau=1 #Variance of center effect on censoring

## DATA SIMULATION
data<-simulate_data(G,Z,prop,beta,theta,cens,pcens,tau)
```

### 📈 Parameters estimation

You can estimate the parameters of the model by using the function
**Parameters_estimation**.

You must specify the method: - Competing risks with shared frailty for
cause 1 (`"CompRisk_frailty"`) - Competing risks without frailty for
cause 1 (`"CompRisk"`) - Standard Cox model with shared frailty
(`"Cox_frailty"`) - Standard Cox model without frailty (`"Cox"`)

You must submit a data frame with at least 3 columns: - times: the
observed time, - status: 0 for right-censoring, i in [1,n] otherwise,
where i is the cause of failure, - clusters: cluster/group indicator. If
there is no cluster, clusters=rep(1,n), with n the size sample.

You also have to specify whether you want to take into account for a
center effect on censoring (cluster_censoring), if you chose
`method = "CompRisk_frailty"`. You may also specify a threshold of
convergence, a threshold for the variance of center effect frailty (in
the case it tends to zero and it might be negligible, hence the
threshold) and a maximum number of iterations.

``` r
method="CompRisk_frailty" #specify the method needed
cluster_censoring=F #center effect on censoring
max_iter=300 #maximum number of iterations
tol=1e-6 #threshold of convergence
threshold=1e-6 #threshold for the variance of center effect frailty

results <- Parameters_estimation(data,method,cluster_censoring,max_iter,tol,threshold)

results$beta #estimation of effects of covariates on failure cause 1
results$theta #(only for method "CompRisk_frailty" and "Cox_frailty") estimation of the variance of center effect frailty
results$u #(only for method "CompRisk_frailty" and "Cox_frailty") estimation of center effect frailty
results$p_value #(only for method "CompRisk_frailty" and "Cox_frailty") p_value for the test theta = 0, if p>0.05 or p is NA, it suggests that the cluster effect may be negligible.
```

### Complexity

The complexity of this algorithm is O(N² + Np + NK + (p + K)³)

where: - N: sample size - T: the number of iterations made, T \<
max_iter, where max_iter is a parameter you can choose - p: the number
of covariables - K: the number of clusters

This package uses sparse matrices with the package `Matrix`, which
improve significantly the speed of the algorithms.

## 🧠 Example on real data

This is under development...

## 📬 Contact

For any questions, issues, or suggestions, please open an issue on
GitHub.

## 🙏 Acknowledgments

Thanks to Lucas Ducrot and Benjamin Delmas for their contribution.

## 🧩 References

The `Matrix` package.

If you use this package in your research, please cite the following
article:

Katsahian S, Resche-Rigon M, Chevret S, Porcher R. Analysing multicentre
competing risks data with a mixed proportional hazards model for the
subdistribution. *Stat Med*. 2006 Dec 30;25(24):4267-78. doi:
[10.1002/sim.2684](https://onlinelibrary.wiley.com/doi/10.1002/sim.2684).
PMID: [16960919.](https://pubmed.ncbi.nlm.nih.gov/16960919/)
