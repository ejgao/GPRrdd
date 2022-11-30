
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
#'
#' @return mean function and covariance function
#' @export
#'
#' @examples
#' X <- matrix(rnorm(24), nrow = 4)
#' y <- c(1, 2, 3, 4)
#' gp_prior(X, y, b = 3, sigma_hat = 0.3, l = 2)
gp_prior <- function(X, Y, b, degree = 0, sigma_hat, choice = 1, l = NULL, alpha = NULL) {
  # y is response, x is matrix of predictors
  # user can specify the degree here
  # first need to split the?plo data by a boundary point/discontinuity point
  # do we need to split by discontinuity point?
  # one dimensional, need to claim X up
  Xc <- X[(1:b - 1), ]
  Yc <- Y[1:b - 1]
  # may have to do regression on one column of X, not every column but need to check on this
  if (degree == 0) {
    mean_function_control <- rep(0, ncol(Xc))
    mean_function_treatment <- rep(0, ncol(Xc))
  } else {
    # mean function for control
    mean_function_control <- mean_function(X, Y, degree)
    # mean function for treatment
    mean_function_treatment <- mean_function(X, Y, degree)
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
  # return means, K
  return(list(mean = mean_function_control, K = K))
}

# mean function
mean_function = function(X, Y, degree){
  mean_func = lm(Y ~ poly(x, degree = degree, raw = TRUE))$coefficients
}

squared_exponential_covfunction <- function(b, X, sigma_hat, l) {
  # b is the discontinuity point
  # for now, allow users to pick sigma_hat, l, will optimize if time permits
  # allow users to choose which function they want to pick
  # allow for calculation of a vector vs a matrix
  n <- length(X)
  # if b isn't supplied, find kernel of X
  if (is.null(b) == TRUE){
    K <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        K[i, j] <- K[j, i] <- sigma_hat * exp({
          -(X[i] - X[j])^2 / (2 * l^2)
        })
      }
    }
  }
  else{
    K <- matrix(0, 1, n)
    for(i in 1:n){
      K[i] = sigma_hat * exp({ -(b - x[i])^2}) / (2*l^2)
    }
  }
  # if b is supplied, find kernel between b and X
  return(K)
}


rational_quad_kernel <- function(X, alpha, l, sigma_hat) {
  n <- length(X)
  if (is.null(b) == TRUE){
    for (i in 1:n) {
      for (j in 1:n) {
        K[i, j] <- K[j, i] <- sigma_hat * (1 + (X[i] - X[j])^2 / (2 * alpha * l^2))^(-alpha)
      }
    }
  }
  else{
    K <- matrix(0, 1, n)
    for(i in 1:n){
      K[i] = sigma_hat * (1 + (b[i] - X[i])^2 / (2 * alpha * l^2))^(-alpha)
    }
  }
  # if b is supplied, find kernel between b and k(1 by n)
  return(K)
}


gp_posterior <- function(X, Y, b, sigma_hat, choice = 1, l = NULL, alpha = NULL) {
  # b is boundary point. it is a scalar, let b the row number you want to let be the boundary
  # mean function at b
  # Xt is stuff after the boundary
  Xc <- X[1:(b - 1)]
  Xt <- X[(b + 1):nrow(n), ]
  Yc <- Y[1:b - 1]
  Yt <- Y[(b + 1):nrow(n)]
  n <- nrow(Xt)
  sigma_y <- var(Y) * diag(n)
  if (choice == 1) {
    Kc <- squared_exponential_covfunction(Xc, sigma_hat, l)
    Kt <- squared_exponential_covfunction(Xt, sigma_hat, l)

  }
  if (choice == 2) {
    Kc <- rational_quad_kernel(Xc, alpha, l, sigma_hat)
    Kt <- rational_quad_kernel(Xt, alpha, l, sigma_hat)
  }
  # need to get posterior mean for both treatment and control
  #diff_t = Yt - mean_function(Xt, Yt, degree)
  #diff_c = Yc - mean_function(Xc, Yc, degree)
  #kernel_t = e + solve(Kt + sigma_y)
  #kernel_c = KTc + solve(Kc + sigma_y)
  # posterior mean is given by mean at the boundary + kernel*(y-mean)
  #post_mean_c = 0 + "quatty" + kernel_t
  #post_mean_treat = 0 + "quatty" + kernel_c
  # next get posterior variance

  # need to get posterior variance for both treatment and control
}



## inherit parameters, add reference
#' Title
#' @inheritParams gp_prior
#' @return
#'
#' @export
#' @examples
create_plot <- function(X, Y, b, choice) {
  # return plot
  # return gp prior plot
  # create prior plot with boundary removed
  Xc <- X[1:(b - 1)]
  Xt <- X[(b + 1):nrow(n)]
  Yc <- Y[1:b - 1]
  Yt <- Y[(b + 1):nrow(n)]
  # call gp_prior, gp_posterior
  abline(v = b)
  prior <- gp_prior(X, Y, b, degree, choice, sigma_hat, alpha)
  # plot fitted values on prior
  # plot fitted values on posterior
  posterior <- gp_posterior(X, Y, b, choice, sigma_hat, alpha)
  # plot
  # return plot
}
