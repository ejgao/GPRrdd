
#' GP Priors
#'
#' @param X matrix input
#' @param Y vector input
#' @param b scalar input should be the row number user wants to remove(has to be greater than 2)
#' @param degree scalar input, default is 0 if no argument given
#' @param sigma_hat vector input
#' @param choice scalar input 1 or 2, default is 1
#' @param l vector input
#' @param alpha vector input
#' @param choice input
#'
#' @return mean function and covariance function
#' @export
#'
#' @examples
#' X <- matrix(rnorm(24), nrow = 4)
#' y <- c(1, 2, 3, 4)
#' gp_prior(X, y, b = 3, sigma_hat = 0.3, l = 2)
gp_prior <- function(X, Y, b, degree = 0, choice = 1, sigma_hat = NULL, l = NULL, alpha = NULL) {
  # y is response, x is matrix of predictors
  # user can specify the degree here
  # first need to split the?plo data by a boundary point/discontinuity point
  Xc <- X[(1:b - 1), ]
  Yc <- Y[1:b - 1]
  # may have to do regression on one column of X, not every column but need to check on this
  if (degree == 0) {
    mean_function_control <- rep(0, ncol(Xc))
    mean_function_treatment <- rep(0, ncol(Xc))
  } else {
    # mean function for control
    mean_function_control <- lm(Yc ~ poly(x = Xc, degree = degree, raw = TRUE))$coefficients
    # mean function for treatment
    mean_function_treatment <- lm(Yc ~ poly(x = Xc, degree = degree, raw = TRUE))$coefficients
  }
  # call covariance function for treatment
  # GP Processes = multivariate normal over a finite set, hence use rmnorm
  # choice of covariance kernel left to user
  if (choice == 1) {
    K <- squared_exponential_covfunction(Xc, sigma_hat, l)
  }
  if (choice == 2) {
    K <- rational_quad_kernel(Xc, alpha, l, sigma_hat)
  }
  # return mean, K
  return(list(mean = mean_function_control, K = K))
}


squared_exponential_covfunction <- function(X, sigma_hat, l) {
  # for now, allow users to pick sigma_hat, l, will optimize if time permits
  # allow users to choose which function they want to pick
  n <- nrow(X)
  K <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- K[j, i] <- sigma_hat * exp({
        -(X[i] - X[j])^2 / (2 * l^2)
      })
    }
  }
  return(K)
}


rational_quad_kernel <- function(X, alpha, l, sigma_hat) {
  n <- nrow(X)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- K[j, i] <- sigma_hat * (1 + (X[i] - X[j])^2 / (2 * alpha * l^2))^(-alpha)
    }
  }
  return(K)
}


gp_posterior <- function(X, Y, b) {
  # b is boundary point. it is a scalar, let b the row number you want to let be the boundary
  # mean function at b
  # Xt is stuff after the boundary
  Xt <- X[(b + 1):nrow(n), ]
  n <- nrow(Xt)
  sigma_y <- var(Y) * diag(n)
  # need to get posterior mean for both treatment and control
  # need to get posterior variance for both treatment and control
}



## inherit parameters, add reference
#' Title
#' @inheritParams gp_prior
#' @return
#'
#' @export
#' @examples
create_plot <- function(b, X, Y) {
  # return plot
  # return gp prior plot
  # create prior plot with boundary removed
  Xc <- X[1:(b - 1)]
  Xt <- X[(b + 1):nrow(n)]
  Yc <- Y[1:b - 1]
  Yt <- Y[(b + 1):nrow(n)]
  # call gp_prior, gp_posterior
  abline(v = X[b])
  prior <- gp_prior(X, Y, b, degree, choice, sigma_hat, alpha)
  # plot fitted values on prior
  # plot fitted values on posterior
  posterior <- gp
  # return plot
}
