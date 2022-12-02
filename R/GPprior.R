
#' GP Priors
#'
#' @param X matrix input
#' @param Y vector input
#' @param b scalar input, should be the discontinuity point in the data
#' @param col_num scalar input allowing user to select particular column of X
#' @param sigma_hat vector input
#' @param degree scalar input, default is NULL if no argument given
#' @param choice scalar input 1 or 2, default is 1
#' @param l scalar input
#' @param alpha vector input
#'
#' @return mean function and covariance function
#' @export
#'
#' @examples
#' X <- matrix(rnorm(24), nrow = 4)
#' y <- c(1, 2, 3, 4)
#' gp_prior(X, y, b = 3, sigma_hat = 0.3, l = 2)
gp_prior <- function(X, Y, b, col_num, sigma_hat, degree = NULL, choice = 1, l = NULL, alpha = NULL) {
  # y is response, x is matrix of predictors
  # user can specify the degree here
  # first need to split the data by a boundary point/discontinuity point
  # do we need to split by discontinuity point?
  # one dimensional, need to claim X up
  X = X[, col_num]
  split_point = which(X == b)
  Xc <- X[1: split_point[1]] # choose first instance data splits
  Yc <- Y[1: split_point[1]]
  Xt <- X[split_point[2]: length(X)]
  Yt <- Y[split_point[2]: length(X)]
  # may have to do regression on one column of X, not every column but need to check on this
  if (is.null(degree) == TRUE) {
    mean_function_control <- rep(0, length(Xc))
    mean_function_treatment <- rep(0, length(Xt))
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
    K <- squared_exponential_covfunction(X, sigma_hat, l)
  }
  if (choice == 2) {
    K <- rational_quad_kernel(X, sigma_hat, l, alpha)
  }
  # return means, K
  return(list(mean_control= mean_function_control, mean_treatment = mean_function_treatment, K = K))
}

# mean function
mean_function = function(X, Y, degree){
  mean_func = lm(Y ~ poly(x, degree = degree, raw = TRUE))$coefficients
}

# squared exponential covariance kernel
squared_exponential_covfunction <- function(X, sigma_hat, l, b = NULL) {
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
    K = matrix(0, 1, n)
    K = sigma_hat * exp({ -(b - X)^2})/ (2 * l^2)
    #K <- matrix(0, 1, n)
    ## vectorize

    # for(i in 1:n){
    #   K[i] = sigma_hat * exp({ -(b - x[i])^2}) / (2*l^2)
    # }
  }
  # if b is supplied, find kernel between b and X
  return(K)
}


# rational quadratic kernel
rational_quad_kernel <- function(X, alpha, l, sigma_hat, b = NULL) {
  n <- length(X)
  if (is.null(b) == TRUE){
    for (i in 1:n) {
      for (j in 1:n) {
        K[i, j] <- K[j, i] <- sigma_hat * (1 + (X[i] - X[j])^2 / (2 * alpha * l^2))^(-alpha)
      }
    }
  }
  else{
    ## vectorize
    K = matrix(0, 1, n)
    K = sigma_hat * (1 + (b - X)^2)/ ((2 * alpha * l^2)^(-alpha))
    # K <- matrix(0, 1, n)
    # for(i in 1:n){
    #   K[i] = sigma_hat * (1 + (b[i] - X[i])^2 / (2 * alpha * l^2))^(-alpha)
    # }
  }
  # if b is supplied, find kernel between b and X(1 by n)
  # cov between X and b is n by 1
  return(K)
}


# gp_posterior
#' Title
#'
#' @param X
#' @param Y
#' @param b
#' @param col_num
#' @param sigma_hat
#' @param choice
#' @param l
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
gp_posterior <- function(X, Y, b, col_num, sigma_hat, degree = NULL, choice = 1, l = NULL, alpha = NULL) {
  # b is boundary point. it is a scalar, let b the row number you want to let be the boundary
  # mean function at b
  # Xt is stuff after the boundary
  X = X[, col_num]
  split_point = which(X == b)
  Xc <- X[1: split_point[1]] # choose first instance data splits
  Yc <- Y[1: split_point[1]]
  Xt <- X[split_point[2]: length(X)]
  Yt <- Y[split_point[2]: length(X)]
  #n <- nrow(Xt)
  # sigma_y <- var(Y) * diag(n)
  if (choice == 1) {
    Kc <- squared_exponential_covfunction(Xc, sigma_hat, l)
    Kt <- squared_exponential_covfunction(Xt, sigma_hat, l)
    Kb_xc <- squared_exponential_covfunction(Xc, sigma_hat, l, b)
    Kb_xt <- squared_exponential_covfunction(Xt, sigma_hat, l, b)
  }
  else if (choice == 2) {
    Kc <- rational_quad_kernel(Xc, alpha, l, sigma_hat)
    Kt <- rational_quad_kernel(Xt, alpha, l, sigma_hat)
    Kb_xc <- rational_quad_kernel(Xc, sigma_hat, l, b)
    Kb_xt <- rational_quad_kernel(Xt, sigma_hat, l, b)
  }
  #need to get posterior mean for both treatment and control
  #diff_t = Yt - mean_function(Xt, Yt, degree)
  #diff_c = Yc - mean_function(Xc, Yc, degree)
  diff_c = Yc
  diff_t = Yt
  kernel_c = Kb_xc %*% solve(Kc + var(y) * diag(length(Yc)))
  kernel_t = Kb_xt %*% solve(Kt + var(y) * diag(length(Yt)))
  #posterior mean is given by mean at the boundary + kernel*(y-mean)
  # case when mean is 0 first
  post_mean_c = 0 + (kernel_c %*% diff_c)
  post_mean_t = 0 + (kernel_t %*% diff_t)
  # next get posterior variance
  posterior_kc = kernel_c %*% t(Kb_xc)
  posterior_kt = kernrel_t %*% t(Kb_xt)
  # need to get posterior variance for both treatment and control
  return(list(posterior_c_mean = post_mean_c, posterior_t_mean = post_mean_t, posterior_c_var = posterior_kc,
              posterior_t_var = posterior_kt))
}


## creating the plot
## inherit parameters, add reference
#' Title
#' @inheritParams gp_prior
#' @return
#'
#' @export
#' @examples
#' @references
#' Branson et al. (2019) A nonparametric Bayesian methodology
#' for regression discontinuity designs,
#' \emph{Journal of Statistical Planning and Inference}
#' \strong{202} 14-30,
#' \doi{10.1016/j.jspi.2019.01.003}
create_plot <- function(X, Y, b, degree, choice, sigma_hat, alpha) {
  # return plot
  # return gp prior plot
  # create prior plot with boundary removed
  Xc <- X[1:(b - 1)]
  Xt <- X[(b + 1):nrow(n)]
  Yc <- Y[1:b - 1]
  Yt <- Y[(b + 1):nrow(n)]
  # call gp_prior, gp_posterior
  prior <- gp_prior(X, Y, b, degree, choice, sigma_hat, alpha)
  # plot fitted values on prior
  # plot fitted values on posterior
  posterior <- gp_posterior(X, Y, b, choice, sigma_hat, alpha)
  # plot
  # return plot
  prior_mean = prior$mean
  post_mean = posterior
  # plot posterior mean
  plot(c(Xc, Xt), c(Yc, Yt))
  abline(v = b)

}
