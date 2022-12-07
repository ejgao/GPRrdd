
#' GP Priors
#'
#' @param Xc vector input for the control group
#' @param Xt vector input for the test group
#' @param Yc vector input for the response of the control group
#' @param Yt vector input for the response of the test group
#' @param sigma_hat scalar input for covariance kernels
#' @param choice scalar input 1 or 2, default is 1, which is the squared covariance kernel. 2 is rational quadratic kernel
#' @param l scalar input for covariance kernel
#' @param alpha scalar input for rational quadratic kernel
#' @param degree scalar input, default is NULL if no argument given
#'
#' @return A list with the elements
#' \item{mean_control}{A nc length vector for the prior mean function for the control group}
#' \item{mean_treatment}{A nt length vector for the prior mean function for the treatment group}
#' \item{cov_c}{An nc by nc matrix for the prior covariance kernel for control group }
#' \item{cov_t}{An nt by nt matrix for the prior covariance kernel for treatment group}
#' @export
#'
#' @examples
#' set.seed(100)
#' xc <- seq(0, 1, length.out = 100)
#' xt <- seq(1, 2, length.out = 100)
#' yc <- 1 * xc + rnorm(100, 0, 0.25)
#' yt <- 2 + 1 * xt + rnorm(100, 0, 0.25)
#' gp_prior(xc, xt, yc, yt, sigma_hat = 1.2, choice = 1, l = 0.7)
#' gp_prior(xc, xt, yc, yt, sigma_hat = 1.2, choice = 2, l = 0.8, alpha = 2)
gp_prior <- function(Xc, Xt, Yc, Yt, sigma_hat, choice = 1, l = NULL, alpha = NULL, degree = NULL) {
  # call mean function for control
  mean_function_control <- mean_function(Xc, Yc, degree)
  # call mean function for treatment
  mean_function_treatment <- mean_function(Xt, Yt, degree)
  # call covariance function for treatment
  # GP Processes = multivariate normal over a finite set
  # choice of covariance kernel left to user
  # 1 or 2 indicates the choice
  if (choice == 1) {
    Kc <- squared_exponential_covfunction(Xc, sigma_hat, l)
    Kt <- squared_exponential_covfunction(Xt, sigma_hat, l)
  }
  else if (choice == 2) {
    Kc <- rational_quad_kernel(Xc, sigma_hat, l, alpha)
    Kt <- rational_quad_kernel(Xt, sigma_hat, l, alpha)
  }
  else{
    stop("Choice should be either 1 or 2")
  }
  # return means, K
  return(list(mean_control = mean_function_control, mean_treatment = mean_function_treatment, cov_c = Kc, cov_t = Kt))
}

# mean function
mean_function <- function(X, Y, degree = NULL) {
  if (is.null(degree) == FALSE) {
    mean_func <- lm(Y ~ poly(X, degree = degree, raw = TRUE))$fitted.values
  } else {
    mean_func <- rep(0, length(X))
  }
}

# squared exponential covariance kernel
squared_exponential_covfunction <- function(X, sigma_hat, l, b = NULL) {
  # b is the discontinuity point
  # for now, allow users to pick sigma_hat, l, will optimize if time permits
  # allow users to choose which function they want to pick
  # allow for calculation of a vector vs a matrix
  n <- length(X)
  # Compatibility check
  if (l <= 0){
    stop("Lengthscale should be positive")
  }
  # if b isn't supplied, find kernel of X
  if (is.null(b) == TRUE) {
    K <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        K[i, j] <- K[j, i] <- sigma_hat * exp({
          -(X[i] - X[j])^2 / (2 * l^2)
        })
      }
    }
  }
  # if b is supplied, find kernel between b and X
  else {
    # vectorize
    K <- as.matrix(sigma_hat * exp({
      -(b - X)^2
    }) / (2 * l^2))
  }
  return(K)
}


# rational quadratic kernel
rational_quad_kernel <- function(X, alpha, l, sigma_hat, b = NULL) {
  n <- length(X)
  # Compatibility Check
  if(l == 0){
    stop("Lengthscale should be positive")
  }
  if(alpha == 0){
    stop("Alpha cannot be zero")
  }
  if (is.null(b) == TRUE) {
    K <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        K[i, j] <- K[j, i] <- sigma_hat * (1 + (X[i] - X[j])^2 / (2 * alpha * l^2))^(-alpha)
      }
    }
  } else {
    ## vectorize
    # if b is supplied, find kernel between b and X(1 by n)
    # cov between X and b is n by 1
    K <- as.matrix((sigma_hat) * (1 + ((b - X)^2 / ((2 * alpha * l^2))))^(-alpha))
  }
  return(K)
}


# gp_posterior
#' Title GP Posterior
#' @inheritParams gp_prior
#'
#' @return A list with the elements
#' \item{posterior_c_mean}{A nc length vector for the posterior mean function for the control group}
#' \item{posterior_t_mean}{A nt length vector for the posterior mean function for the treatment group}
#' \item{posterior_c_var}{An nc by nc matrix for the posterior covariance kernel for control group}
#' \item{posterior_t_var}{An nc by nc matrix for the posterior covariance kernel for treatment group}
#'
#' @export
#' @examples
#' set.seed(100)
#' xc <- seq(0, 1, length.out = 100)
#' xt <- seq(1, 2, length.out = 100)
#' yc <- 1 * xc + rnorm(100, 0, 0.25)
#' yt <- 2 + 1 * xt + rnorm(100, 0, 0.25)
#' gp_posterior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1.2, choice = 1, l = 0.7)
gp_posterior <- function(Xc, Xt, Yc, Yt, sigma_hat, choice = 1, l = NULL, alpha = NULL, degree = NULL) {
  # mean function at b
  # store Kb_xc, Kb_xt as vectors
  post_mean_c <- rep(NA, length(Xc))
  post_mean_t <- rep(NA, length(Xt))
  var_c <- rep(NA, length(Xc))
  var_t <- rep(NA, length(Xt))
  ##########################
  # y - mean of x for both control and treatment groups
  # get means for Xt and Xc
  meansb_c <- mean_function(Xc, Yc, degree)
  meansb_t <- mean_function(Xt, Yt, degree)
  diff_t <- Yt - mean_function(Xt, Yt, degree)
  diff_c <- Yc - mean_function(Xc, Yc, degree)
  # two choices for cov kernel
  # denote var of y as var(c(Yc, Yt))
  if (choice == 1) {
    Kc <- squared_exponential_covfunction(Xc, sigma_hat, l)
    Kt <- squared_exponential_covfunction(Xt, sigma_hat, l)
    kernel_c <- solve(Kc + var(c(Yc, Yt)) * diag(length(Yc)))
    kernel_t <- solve(Kt + var(c(Yc, Yt)) * diag(length(Yt)))
    # looping around all the control groups
    for (i in 1:length(Xc)) {
      post_mean_c[i] <- meansb_c[i] + t(squared_exponential_covfunction(Xc, sigma_hat, l, Xc[i])) %*% kernel_c %*% diff_c
      var_c[i] <- t(squared_exponential_covfunction(Xt, sigma_hat, l, Xc[i])) %*% kernel_c %*% squared_exponential_covfunction(Xt, sigma_hat, l, Xc[i])
    }
    # looping around all the treatment groups
    for (i in 1:length(Xt)) {
      post_mean_t[i] <- meansb_t[i] + t(squared_exponential_covfunction(Xt, sigma_hat, l, Xt[i])) %*% kernel_t %*% diff_t
      var_t[i] <- t(squared_exponential_covfunction(Xt, sigma_hat, l, Xt[i])) %*% kernel_t %*% squared_exponential_covfunction(Xt, sigma_hat, l, Xt[i])
    }
  } else if (choice == 2) {
    Kc <- rational_quad_kernel(Xc, alpha, l, sigma_hat)
    Kt <- rational_quad_kernel(Xt, alpha, l, sigma_hat)
    kernel_c <- solve(Kc + var(c(Yc, Yt)) * diag(length(Yc)))
    kernel_t <- solve(Kt + var(c(Yc, Yt)) * diag(length(Yt)))
    # looping around all the control groups
    for (i in 1:length(Xc)) {
      post_mean_c[i] <- meansb_c[i] + t(rational_quad_kernel(Xc, alpha, l, sigma_hat, Xc[i])) %*% kernel_c %*% diff_c
      var_c[i] <- t(rational_quad_kernel(Xc, alpha, l, sigma_hat, Xc[i])) %*% kernel_c %*% rational_quad_kernel(Xc, alpha, l, sigma_hat, Xc[i])
    }
    # looping around all the treatment groups
    for (i in 1:length(Xt)) {
      post_mean_t[i] <- meansb_t[i] + t(rational_quad_kernel(Xt, alpha, l, sigma_hat, Xt[i])) %*% kernel_t %*% diff_t
      var_t[i] <- t(rational_quad_kernel(Xt, alpha, l, sigma_hat, Xt[i])) %*% kernel_t %*% rational_quad_kernel(Xt, alpha, l, sigma_hat, Xt[i])
    }
  }
  else{
    stop("Choice should be either 1 or 2")
  }
  return(list(
    posterior_c_mean = post_mean_c, posterior_t_mean = post_mean_t, posterior_c_var = var_c,
    posterior_t_var = var_t
  ))
}


## creating the plot
## inherit parameters, add reference
#' Creating RDD Plot
#' @param X matrix input
#' @param Y matrix input
#' @param b discontinuity point
#' @param col_num column of X that the user wants in order to do regression discontinuty design on
#' @param sigma_gp scalar input--indicates the variance of the Gaussian Process model
#' @inheritParams gp_prior
#' @return RDD Plot as well as estimated treatment effect
#'
#' @export
#' @examples
#' set.seed(100)
#' sc <- seq(0, 1, length.out = 100)
#' st <- seq(1, 2, length.out = 100)
#' yc <- 1 * sc + rnorm(100, 0, 0.25)
#' yt <- 2 + 1 * st + rnorm(100, 0, 0.25)
#' x <- c(sc, st)
#' x <- as.matrix(x)
#' y <- c(yc, yt)
#' create_plot(X = x, Y = y, b = 1, col_num = 1, sigma_gp = 2, sigma_hat = 1.2, choice = 1, l = 0.7)
#' @references
#' Branson et al. (2019) A nonparametric Bayesian methodology
#' for regression discontinuity designs,
#' \emph{Journal of Statistical Planning and Inference}
#' \strong{202} 14-30,
#' \doi{10.1016/j.jspi.2019.01.003}
create_plot <- function(X, Y, b, col_num, sigma_gp, sigma_hat, choice = 1, l = NULL, alpha = NULL, degree = NULL) {
  # return plot
  # return gp prior plot
  # create prior plot with boundary removed
  X <- X[, col_num]
  # split at discontinuity point
  split_point <- which(X == b)
  # choose first instance where there is a data split
  Xc <- X[1:split_point[1]]
  Yc <- Y[1:split_point[1]]
  # split y also
  Xt <- X[split_point[2]:length(X)]
  Yt <- Y[split_point[2]:length(X)]
  # call gp_prior, gp_posterior
  prior <- gp_prior(Xc, Xt, Yc, Yt, sigma_hat, choice, l, alpha, degree)
  # plot fitted values on prior
  # plot fitted values on posterior
  posterior <- gp_posterior(Xc, Xt, Yc, Yt, sigma_hat, choice, l, alpha, degree)
  # obtain posterior means
  post_mean_c <- posterior$posterior_c_mean
  post_mean_t <- posterior$posterior_t_mean
  # obtain posterior covariances
  # need sigma_gp
  sigma_var_c <- posterior$posterior_c_var
  sigma_var_t <- posterior$posterior_t_var
  # ideally, would compute sigma_gp through MLE
  # However, if we let user choose, we have to ensure that the difference between sigma_gp and posterior sigma is at least 0
  if (sigma_gp < max(sigma_var_c) | sigma_gp < max(sigma_var_t)) {
    sigma_gp <- max(max(sigma_var_c), max(sigma_var_t))
  }
  posterior_var_c <- sigma_gp - sigma_var_c
  posterior_var_t <- sigma_gp - sigma_var_t
  # plot original data
  plot(c(Xc, Xt), c(Yc, Yt), main = "RDD Plot", xlab = "x", ylab = "y")
  # graph discontinuity point by using vertical line through that point
  abline(v = b)
  # graph the fit
  lines(Xc, post_mean_c, col = "red")
  lines(Xt, post_mean_t, col = "blue")
  # plot 95% confidence bands
  # confidence bands for control groups
  lower_c <- post_mean_c - (qt(0.95, length(Xc) - 1) * sqrt(posterior_var_c)) / (sqrt(length(Xc)))
  upper_c <- post_mean_c + (qt(0.95, length(Xc) - 1) * sqrt(posterior_var_c)) / sqrt(length(Xc))
  # confidence bands for treatment groups
  lower_t <- post_mean_t - (qt(0.95, length(Xt) - 1) * sqrt(posterior_var_t)) / sqrt(length(Xt))
  upper_t <- post_mean_t + (qt(0.95, length(Xt) - 1) * sqrt(posterior_var_t)) / sqrt(length(Xt))
  # graphing bounds for Xc
  lines(Xc, lower_c, col = "red", lty = 5)
  lines(Xc, upper_c, col = "red", lty = 5)
  # graphing bounds for Xt
  lines(Xt, lower_t, col = "blue", lty = 5)
  lines(Xt, upper_t, col = "blue", lty = 5)
  # add a legend to differentiate between every line
  # find treatment effect(scalar)
  post_treatment_effect <- post_mean_t[1] - post_mean_c[length(Xc)]
  # Calculate the 95% CI around post_treatment effect
  ci_effect_upper <- post_treatment_effect + 1.96 * sqrt((posterior_var_t[1] + posterior_var_c[length(Xc)]))
  ci_effect_lower <- post_treatment_effect - 1.96 * sqrt((posterior_var_t[1] + posterior_var_c[length(Xc)]))
  ci = c(ci_effect_lower, ci_effect_upper)
  # return treatment effect(scalar)
  return(list(post_treatment_effect = post_treatment_effect, ci = ci))
}
