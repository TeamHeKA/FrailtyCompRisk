# multicenterCompRisk

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/multicenterCompRisk)](https://CRAN.R-project.org/package=multicenterCompRisk)
[![R-CMD-check](https://github.com/DELMASben/multicenterCompRisk/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DELMASben/multicenterCompRisk/actions/workflows/R-CMD-check.yaml)
[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

The goal of **`multicenterCompRisk`** is to provide tools for **competing risks survival analysis** in **multicenter studies**, accounting for **unobserved heterogeneity** (frailty) at the center level as proposed in the following article:

Katsahian S, Resche-Rigon M, Chevret S, Porcher R. Analysing multicentre competing risks data with a mixed proportional hazards model for the subdistribution. *Stat Med*. 2006 Dec 30;25(24):4267-78. doi: [`10.1002/sim.2684`]{https://onlinelibrary.wiley.com/doi/10.1002/sim.2684}. PMID: [`16960919.`]{https://pubmed.ncbi.nlm.nih.gov/16960919/}

This package is particularly useful for statisticians and epidemiologists analyzing time-to-event data where individuals are nested within centers (e.g., hospitals or clinics), and where multiple causes of failure may occur.

It includes functions to:
- Simulate clustered competing risks data,
- Fit cause-specific Cox models with shared frailty, following the method presented by Sandrine Katsahian in "Analysing multicentre competing risks data with a mixed proportional hazards model for the subdistribution."
- Fit cause-specific Cox models,
- Fit Cox models with shared frailty,
- Fit Cox models.

The REML method is used when shared frailty is modeled; otherwise, the ML method is applied.

## 📦 Installation

You can install the development version of `multicenterCompRisk` from GitHub using either [`pak`](https://pak.r-lib.org/) or [`devtools`](https://github.com/r-lib/devtools):

```r
# Using pak (recommended)
# install.packages("pak")
pak::pak("DELMASben/multicenterCompRisk")

# Or using devtools
# install.packages("devtools")
devtools::install_github("DELMASben/multicenterCompRisk")
```

