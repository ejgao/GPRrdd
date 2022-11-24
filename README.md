
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Gaussian Processes Regression on Regression Discontinuity Designs Package

The following is a short description explaning the use of the package and how the package can be installed.

<!-- badges: start -->
<!-- badges: end -->

## Intended Use

In a regression discontinuity design (RDD), we have a discontinuity
point at a particular covariate value, for which after that point, the
treatment is assigned. In other words, we have a control (before the
discontinuity point) and a treatment (everything after the discontinuity
point). The main goal of the RDD is to estimate the treatment effect.
This is usually done through local linear regression (LLR), where the
main goal is to estimate the treatment effect at the boundary (the
discontinuity point). However, it has been discovered that LLR can have
poor inferential properties. This motivates the non-parametric approach
used in this paper \[Branson et al., 2019\], which employs Gaussian
process regression (GPR) to find an initial prior on the mean response
functions for both our control and treatment groups. After, inference is
done to estimate the treatment effect. This package will show how we can
use Gaussian processes on RDDs.

## Installation

You can install the development version of GPRdd from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("ejgao/GPRrdd")
```

Upon installing, please run

``` r
library(GPRrdd)
```

in order to gain access to the three functions that are available:
gp_prior, gp_posterior, and create_plot.

## What is Left

The actual coding portion of the package is nearly complete and I plan
on finishing up create_plot and gp_posterior. I then plan on creating a
vignette with examples. Furthermore, I also need to document some of my
functions and also polish up all my functions. I know my functions work;
however, I am not sure they work properly so I will be checking on that,
