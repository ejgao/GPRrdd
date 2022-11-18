
#' GP Priors
#'
#' @param X matrix input
#' @param Y vector input
#' @param b scalar input should be the row number user wants to remove
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
gp_prior = function(X, Y, b, degree = 0, choice = 1, sigma_hat = NULL, l = NULL, alpha = NULL){
  # y is response, x is matrix of predictors
  # user can specify the degree here
  # first need to split the data by a boundary point/discontinuity point
  Xc = X[(1:b-1), ]
  if (degree == 0){
    mean_function_control = rep(0, ncol(Xc))
    mean_function_treatment = rep(0, ncol(Xc))
  }
  else{
    # mean function for control
    mean_function_control <- lm(Y ~ poly(Xc,degree))$coefficients
    # mean function for treatment
    mean_function_treatment <- lm(Y ~ poly(Xc,degree))$coefficients
  }
  # call covariance function for treatment
  # GP Processes = multivariate normal over a finite set, hence use rmnorm
  # choice of covariance kernel left to user
    if (choice == 1){
      K = squared_exponential_covfunction(Xc, sigma_hat, l)
    }
    if (choice == 2){
      K = rational_quad_kernel(Xc, alpha, l, sigma_hat)
    }
  # return mean, K
  return(list(mean = mean_function_control, K = K))
}


#' Title
#'
#' @param X
#' @param sigma_hat
#' @param l
#'
#' @return
#' @export
#'
#' @examples
squared_exponential_covfunction = function(X, sigma_hat, l){
  # for now, allow users to pick sigma_hat, l, will optimize if time permits
  # allow users to choose which function they want to pick
  n = nrow(X)
  K = matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] = K[j, i] = sigma_hat * exp({-(X[i] - X[j])^2/(2 * l^2)})
    }
  }
  return(K)
}

rational_quad_kernel = function(X, alpha, l, sigma_hat){
  n = nrow(X)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] = K[j, i] = sigma_hat * (1 + (X[i] - X[j])^2 / (2 * alpha * l^2))^(-alpha)
    }
  }
  return(K)
}


gp_posterior = function(X, Y, b){
  # b is boundary point. it is a scalar, let b the row number you want to let be the boundary
  # mean function at b
  # Xt is stuff after the boundary
  Xt = X[(b + 1): nrow(n), ]
  n = nrow(Xt)
  sigma_y = var(Y) * diag(n)
}

# need a readMe

create_plot = function(b, X, Y){
  # return plot
  # return gp prior plot
  # create prior plot with boundary removed

  return(1)
}






