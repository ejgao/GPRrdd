
#' GP Priors
#'
#' @param degree scalar input, default is 2 if no argument given
#' @param X matrix input
#' @param y vector input
#'
#' @return multivariate normal distribution with given mean function and covariance function
#' @export
#'
#' @examples
#' gp_prior(degree = 3, X, y)
gp_prior = function(degree = 0, X, y){
  # y is response, x is matrix of predictors
  # user can specify the degree here
  if (degree == 0){
    mean_function_control = rep(0, ncol(X))
    mean_function_treatment = rep(0, ncol(X))
  }
  else{
    mean_function_control <- lm(Y ~ poly(X,degree))
    # mean function for treatment
    mean_function_treatment <- lm(Y ~ poly(X,degree))
  }
  # first find the distance between every point in x, then take exponential of that
  # call covariance function for treatment
  cov = covariance_function(X,y, signa_hat, l)
  K = cov$cov_matrix
  # allow user to specify the number of observations
  # GP Processes = multivariate normal over a finite set, hence use rmnorm
  prior_control_group = rmvnorm(mean_function_control, K) # fix these
  prior_treatment_group = rmvnorm(mean_function_treatment, K)
  # Apply rmvnorm on both control and treatment
  # return list
  return(2)
}

covariance_function = function(X, y, sigma_hat, l){
  # for now, allow users to pick sigma_hat, l, will optimize if time permits
  # allow users to choose which function they want to pick
  # compute distance between any two rows of X
  # randomize to obtain two rows of X
  # this covariance function does not change for control and treatment groups.
  # first covariance function is squared exponential
  random_indices = sample(nrow(X), 2, replace = F)
  distance_matrix = X[random_indices, ]
  x1 = distance_matrix[1,]
  x2 = distance_matrix[2,]
  # take the distance between two matrices
  # let user choose covariance matrix
  # list of options
  # https://towardsdatascience.com/gaussian-process-kernels-96bafb4dd63e#:~:text=Perhaps%20the%20most%20widely%20used,L%20the%20kernel%20length%20scale.
  total_d = crossprod(X1, X2)
  cov_matrix = sigma_hat * exp(total_d/l)
  return(cov_matrix)
}

gp_posterior = function(X, y){

}



create_plot = function(b, X, y){
  # return plot
  return(1)
}

