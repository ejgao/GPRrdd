
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
  if (choice == 2) {
    Kc <- rational_quad_kernel(Xc, sigma_hat, l, alpha)
    Kt <- rational_quad_kernel(Xt, sigma_hat, l, alpha)

  }
  # return means, K
  return(list(mean_control= mean_function_control, mean_treatment = mean_function_treatment, cov_c = Kc, cov_t = Kt))
}

# mean function
mean_function = function(X, Y, degree = NULL){
  if (is.null(degree) == FALSE){
    mean_func = lm(Y ~ poly(x, degree = degree, raw = TRUE))$fitted.values
  }
  else{
    mean_func = rep(0, length(X))
  }
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
  # if b is supplied, find kernel between b and X
  else{
    # vectorize
    K = as.matrix(sigma_hat * exp({ -(b - X)^2})/ (2 * l^2))
  }
  return(K)
}


# rational quadratic kernel
rational_quad_kernel <- function(X, alpha, l, sigma_hat, b = NULL) {
  n <- length(X)
  if (is.null(b) == TRUE){
    K <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        K[i, j] <- K[j, i] <- sigma_hat * (1 + (X[i] - X[j])^2 / (2 * alpha * l^2))^(-alpha)
      }
    }
  }
  else{
    ## vectorize
    # if b is supplied, find kernel between b and X(1 by n)
    # cov between X and b is n by 1
    K = as.matrix((sigma_hat) * (1 + ((b - X)^2 / ((2 * alpha * l^2))))^(-alpha))
  }
  return(K)
}


# gp_posterior
#' Title
#' @inheritParams gp_prior
#'
#' @return
#' @export
#'
#' @examples
gp_posterior <- function(Xc, Xt, Yc, Yt, sigma_hat, choice = 1, l = NULL, alpha = NULL, degree = NULL) {
  # mean function at b
  # store Kb_xc, Kb_xt as vectors
  post_mean_c = rep(NA, length(Xc))
  post_mean_t = rep(NA, length(Xt))
  var_c = rep(NA, length(Xc))
  var_t = rep(NA, length(Xt))
  ##########################
  # y - mean of x for both control and treatment groups
  diff_t = Yt - mean_function(Xt, Yt, degree)
  diff_c = Yc - mean_function(Xc, Yc, degree)
  # two choices for cov kernel
  # denote var of y as var(c(Yc, Yt))
  if (choice == 1) {
    Kc <- squared_exponential_covfunction(Xc, sigma_hat, l)
    Kt <- squared_exponential_covfunction(Xt, sigma_hat, l)
    kernel_c = solve(Kc + var(c(Yc, Yt)) * diag(length(Yc)))
    kernel_t = solve(Kt + var(c(Yc, Yt)) * diag(length(Yt)))
    # looping around all the control groups
    for(i in 1:length(Xc)){
      post_mean_c[i] <- t(squared_exponential_covfunction(Xc, sigma_hat, l, Xc[i])) %*% kernel_c %*% diff_c
      var_c[i] <- t(squared_exponential_covfunction(Xt, sigma_hat, l, Xc[i])) %*% kernel_c %*% squared_exponential_covfunction(Xt, sigma_hat, l, Xc[i])
    }
    # looping around all the treatment groups
    for(i in 1:length(Xt)){
      post_mean_t[i] <- t(squared_exponential_covfunction(Xt, sigma_hat, l, Xt[i])) %*% kernel_t %*% diff_t
      var_t[i] <- t(squared_exponential_covfunction(Xt, sigma_hat, l, Xt[i])) %*% kernel_t %*% squared_exponential_covfunction(Xt, sigma_hat, l, Xt[i])
      }
  }
  else if (choice == 2) {
    Kc <- rational_quad_kernel(Xc, alpha, l, sigma_hat)
    Kt <- rational_quad_kernel(Xt, alpha, l, sigma_hat)
    kernel_c = solve(Kc + var(c(Yc, Yt)) * diag(length(Yc)))
    kernel_t = solve(Kt + var(c(Yc, Yt)) * diag(length(Yt)))
    # looping around all the control groups
    for(i in 1:length(Xc)){
      post_mean_c[i] <- t(rational_quad_kernel(Xc, alpha, l, sigma_hat, Xc[i])) %*% kernel_c %*% diff_c
      var_c[i] = t(rational_quad_kernel(Xc, alpha, l, sigma_hat, Xc[i])) %*% kernel_c %*% rational_quad_kernel(Xc, alpha, l, sigma_hat, Xc[i])

    }
    # looping around all the treatment groups
    for(i in 1:length(Xt)){
      post_mean_t[i] <- t(rational_quad_kernel(Xt, alpha, l, sigma_hat, Xt[i])) %*% kernel_t %*% diff_t
      var_t[i] = t(rational_quad_kernel(Xt, alpha, l, sigma_hat, Xt[i])) %*% kernel_t %*% rational_quad_kernel(Xt, alpha, l, sigma_hat, Xt[i])
    }
  }
  return(list(posterior_c_mean = post_mean_c, posterior_t_mean = post_mean_t, posterior_c_var = var_c,
              posterior_t_var = var_t))
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
create_plot <- function(X, Y, b, col_num, sigma_gp, sigma_hat, choice = 1, l = NULL, alpha = NULL, degree = NULL) {
  # return plot
  # return gp prior plot
  # create prior plot with boundary removed
  X = X[, col_num]
  # split at discontinuity point
  split_point = which(X == b)
  # choose first instance where there is a data split
  Xc <- X[1: split_point[1]]
  Yc <- Y[1: split_point[1]]
  # split y also
  Xt <- X[split_point[2]: length(X)]
  Yt <- Y[split_point[2]: length(X)]
  # call gp_prior, gp_posterior
  prior <- gp_prior(Xc, Xt, Yc, Yt, sigma_hat, choice, l, alpha, degree)
  # plot fitted values on prior
  # plot fitted values on posterior
  posterior <- gp_posterior(Xc, Xt, Yc, Yt, sigma_hat, choice, l, alpha, degree)
  # obtain posterior means
  post_mean_c = posterior$posterior_c_mean
  post_mean_t = posterior$posterior_t_mean
  # obtain posterior covariances
  # need sigma_gp
  sigma_var_c = posterior$posterior_c_var
  sigma_var_t = posterior$posterior_t_var
  # ideally, would compute sigma_gp through MLE
  # However, if we let user choose, we have to ensure that the difference between sigma_gp and posterior sigma is at least 0
  if(sigma_gp < max(sigma_var_c) | sigma_gp < max(sigma_var_t)){
    sigma_gp = max(max(sigma_var_c), max(sigma_var_t))
  }
  posterior_var_c = sigma_gp - sigma_var_c
  posterior_var_t = sigma_gp - sigma_var_t
  # plot original data
  plot(c(Xc, Xt), c(Yc, Yt), main = "RDD Plot", xlab = "x", ylab = "y")
  # graph discontinuity point by using vertical line through that point
  abline(v = b)
  # graph the fit
  lines(Xc, post_mean_c, col = "red")
  lines(Xt, post_mean_t, col = "blue")
  # plot 95% confidence bands
  # confidence bands for control groups
  lower_c = post_mean_c - (qt(0.95, length(Xc) - 1) * posterior_var_c)/(sqrt(length(Xc)))
  upper_c = post_mean_c + (qt(0.95, length(Xc) - 1) * posterior_var_c)/sqrt(length(Xc))
  # confidence bands for treatment groups
  lower_t = post_mean_t - (qt(0.95, length(Xt) - 1) * posterior_var_t)/sqrt(length(Xt))
  upper_t = post_mean_t + (qt(0.95, length(Xt) - 1) * posterior_var_t)/sqrt(length(Xt))
  # graphing bounds for Xc
  lines(Xc, lower_c, col = "red", lty = 5)
  lines(Xc, upper_c, col = "red", lty = 5)
  #graphing bounds for Xt
  lines(Xt, lower_t, col = "blue", lty = 5)
  lines(Xt, upper_t, col = "blue", lty = 5)
  # add a legend to differentiate between every line
  # legend("topright",inset = c(-.2, 0.2), legend = c("Fit", "95% CI", "Fit", "95% CI"), col = c("red", "red", "blue", "blue"),
  #        lty = c(1,5,1,5))
  # find treatment effect(scalar)
  post_treatment_effect = post_mean_t[1] - post_mean_c[length(Xc)]
  # return treatment effect(scalar)
  return(post_treatment_effect)

}
