
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GPRrdd

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/ejgao/GPRrdd/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ejgao/GPRrdd?branch=master)
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
# install with vignette
devtools::install_github("ejgao/GPRrdd", build_vignettes = TRUE)
```

Upon installing, please run

``` r
library(GPRrdd)
```

in order to gain access to the three functions that are available:
gp_prior, gp_posterior, and create_plot.

## Small Example

Now, an example of create_plot will be shown. Since gp_prior and
gp_posterior have long outputs, their usage be shown; however, the
output will not be displayed. More information can be found in the
vignette, however.

``` r
set.seed(100)
sc = seq(0, 1, length.out = 50)
st = seq(1, 2, length.out = 50)
yc = 1 * sc + rnorm(50, 0, 0.25)
yt = 2 + 1*st + rnorm(50, 0, 0.25)
x = c(sc, st)
x = as.matrix(x)
y = c(yc, yt)
```

Then, to use gp_prior, we can do:

``` r
gp_prior(Xc = sc, Xt = st, Yc = yc, Yt = yt, sigma_hat = 1.2, l = 0.7)
```

Similarly, to use gp_posterior, we can do:

``` r
gp_posterior(Xc = sc, Xt = st, Yc = yc, Yt = yt, sigma_hat = 1.2, l = 0.7)
```

Finally, we can use create_plot to return the model fit on both the
control and treatment groups, the estimated treatment effects, as well
as the confidence interval around the estimated treatment effects.

``` r
create_plot(X=x, Y=y, b = 1, col_num = 1, sigma_gp = 2, sigma_hat = 1.2, choice = 1, l = 0.7)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

    #> $post_treatment_effect
    #> [1] 1.886551
    #> 
    #> $ci
    #> [1] -0.7697033  4.5428051

## References

Branson et al.(2019). A nonparametric Bayesian methodology for
regression discontinuity designs. Journal of Statistical Planning and
Inference, 202, 14-30. <https://doi.org/10.1016%2Fj.jspi.2019.01.003>
