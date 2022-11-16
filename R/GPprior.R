
#' GP Priors
#'
#' @param degree scalar input, default is 2 if no argument given
#' @param X matrix input
#' @param y vector input
#' @param choice input, choose 1 if
#'
#' @return multivariate normal distribution with given mean function and covariance function
#' @export
#'
#' @examples
#' gp_prior(degree = 3, X, y)
#' gp_prior(X, y)
gp_prior = function(degree = 0, X, Y, b, choice){
  # y is response, x is matrix of predictors
  # user can specify the degree here
  # first need to split the data by a boundary point.
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
  K = covariance_function(choice)
  # allow user to specify the number of observations
  # GP Processes = multivariate normal over a finite set, hence use rmnorm
  # prior_control_group = rmvnorm(mean_function_control, K) # fix these
  # prior_treatment_group = rmvnorm(mean_function_treatment, K)
  # return list
  return(list(prior_control = prior_control_group, prior_treatment = prior_treatment_group))
}


squared_exponential_covfunction = function(X, sigma_hat, l){
  # for now, allow users to pick sigma_hat, l, will optimize if time permits
  # allow users to choose which function they want to pick
  n = nrow(x)
  K = matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] = K[j, i] = sigma_hat * exp({-(x[i] - x[j])^2/(2 * l^2)})
    }
  }
  return(K)
}

rational_quad_kernel = function(X, alpha, l, sigma_hat){
  n = nrow(x)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] = K[j, i] = sigma_hat * (1 + (x[i] - x[j])^2 / (2 * alpha * l^2))^(-alpha)
    }
  }
  return(K)
}

covariance_function = function(choice){
  if (choice == 1){
    K = squared_exponential_covfunction(X, Y, sigma_hat, l)
  }
  if (choice == 2){
    K = rational_quad_kernel(X, alpha, l, sigma_hat)
  }
  return(K)
}

gp_posterior = function(X, Y, b){
  # b is boundary point. it is a scalar, let b the row number you want to let be the boundary
  # mean function at b
  # Xt is stuff after the boundary
  n = length(n)
  Xt = X[(b + 1): nrow(n), ]
  posterior_control = cov(Xt, Xt)
  sigma_y = var(y) * diag(n)
}

# need a readMe

create_plot = function(b, X, Y){
  # return plot
  return(1)
}


