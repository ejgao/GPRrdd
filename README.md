
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GPRdd

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
# install without vignette
devtools::install_github("ejgao/GPRrdd")
#> Skipping install of 'GPRdd' from a github remote, the SHA1 (f6e6d5d8) has not changed since last install.
#>   Use `force = TRUE` to force installation
# install with vignette
devtools::install_github("ejgao/GPRrdd", build_vignettes = TRUE)
#> Skipping install of 'GPRdd' from a github remote, the SHA1 (f6e6d5d8) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

Upon installing, please run

``` r
library(GPRdd)
```

in order to gain access to the three functions that are available:
gp_prior, gp_posterior, and create_plot.

## Small Example

Although all three functions are available, the main function is
create_plot. A short example is illustrated below:

``` r
set.seed(100)
sc = seq(0, 1, length.out = 100)
st = seq(1, 2, length.out = 100)
yc = 1 * sc + rnorm(100, 0, 0.25)
yt = 2 + 1*st + rnorm(100, 0, 0.25)
x = c(sc, st)
x = as.matrix(x)
y = c(yc, yt)
create_plot(X=x, Y=y, b = 1, col_num = 1, sigma_gp = 2, sigma_hat = 1.2, choice = 1, l = 0.7)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

    #> [1] 1.93442

## What is Left

The actual coding portion of the package is nearly complete and I plan
on finishing up create_plot and gp_posterior. I then plan on creating a
vignette with examples. Regarding the functions for the package, they
are almost all done. The functions gp_posterior and create_plot need
some work. Furthermore, I also need to document some of my functions and
also polish up all my functions. I know the functions that are complete
work; however, I am not sure they work properly so I will be checking on
that.

## References

Branson et al., 2019\] Branson, Z., Rischard, M., Bornn, L., & Miratrix,
W, L (2019). A nonparametric Bayesian methodology for regression
discontinuity designs. Journal of Statistical Planning and Inference,
202, 14-30. <https://doi.org/10.1016%2Fj.jspi.2019.01.003>
