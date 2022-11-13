
#' GP Priors
#'
#' @param degree
#' @param covariance
#' @param X
#' @param y
#' @param n
#'
#' @return
#' @export
#'
#' @examples
gp_prior = function(degree = 2, covariance, X, y, n){
  # y is response, x is matrix of predictors
  # user can specify the degree here
  mean_function_control <- lm(Y ~ poly(X,degree))
  # mean function for treatment
  mean_function_treatment <- lm(Y ~ poly(X,degree))
  # first find the distance between every point in x, then take exponential of that
  # call covariance function for treatment
  covariance_function(X,y, signa_hat, l)
  # allow user to specify the number of observations
  # GP Processes = multivariate normal over a finite set, hence use rmnorm
  #rmvnorm()
  # Apply rmvnorm on both control and treatment
  return(1)
}

#' Title
#'
#' @param X
#' @param y
#' @param sigma_hat
#' @param l
#'
#' @return
#' @export
#'
#' @examples
covariance_function = function(X, y, sigma_hat, l){
  # for now, allow users to pick sigma_hat, l, will optimize if time permits
  # compute distance between any two rows of X
  # randomize to obtain two rows of X
  # this covariance function does not change for control and treatment groups.
  random_indices = sample(nrow(X), 2, replace = F)
  distance_matrix = X[random_indices, ]
  x1 = distance_matrix[1,]
  x2 = distance_matrix[2,]
  # take the distance between two matrices
  total_d = crossprod(X1, X2)
  cov_matrix = sigma_hat * exp(total_d/l)
}




create_plot = function(b, X, y){
  # return plot
}
