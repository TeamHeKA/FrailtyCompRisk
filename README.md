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
library(frailtycomprisk)

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

## 🧠 Example on real data

This is under development...

## 📬 Contact

For any questions, issues, or suggestions, please open an issue on
GitHub.

## 🙏 Acknowledgments

Thanks to Lucas Ducrot and Benjamin Delmas for their contribution.

## 🧩 References

If you use this package in your research, please cite the following
article:

Katsahian S, Resche-Rigon M, Chevret S, Porcher R. Analysing multicentre
competing risks data with a mixed proportional hazards model for the
subdistribution. *Stat Med*. 2006 Dec 30;25(24):4267-78. doi:
[10.1002/sim.2684](https://onlinelibrary.wiley.com/doi/10.1002/sim.2684).
PMID: [16960919.](https://pubmed.ncbi.nlm.nih.gov/16960919/)
