
## Test that dimensions in gp_prior are correct
set.seed(100)
xc <- seq(0, 1, length.out = 5)
xt <- seq(1, 2, length.out = 5)
yc <- 1 * xc + rnorm(5, 0, 0.25)
yt <- 2 + 1 * xt + rnorm(5, 0, 0.25)
x <- as.matrix(c(xc, xt))
y <- c(yc, yt)

test_that("dimensions in each element of gp_prior are correct", {
  expect_equal(length(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$mean_control), length(xc))
  expect_equal(length(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1, degree = 1)$mean_control), length(xc))
  expect_equal(length(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$mean_treatment), length(xt))
  expect_equal(dim(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$cov_c)[1], length(xc))
  expect_equal(dim(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$cov_c)[2], length(xc))
  expect_equal(dim(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$cov_t)[1], length(xt))
  expect_equal(dim(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$cov_t)[2], length(xt))
  expect_equal(dim(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 2, l = 1, alpha = 2)$cov_c)[1], length(xc))

})

test_that("dimensions in each element of gp_posterior are correct", {
  expect_equal(length(gp_posterior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$posterior_c_mean), length(xc))
  expect_equal(length(gp_posterior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$posterior_t_mean), length(xt))
  expect_equal(length(gp_posterior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$posterior_c_var), length(xc))
  expect_equal(length(gp_posterior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$posterior_t_var), length(xc))
  expect_equal(length(gp_posterior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 2, l = 1, alpha = 1)$posterior_t_var), length(xc))
  expect_equal(length(gp_posterior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 2, l = 2, alpha = 1.3)$posterior_t_var), length(xc))
})


test_that("dimensions in each element of create_plot are correct", {
  expect_equal(length(create_plot(X = x, Y = y, b = 1, col_num = 1, sigma_gp = 2, sigma_hat = 1.2, choice = 1, l = 0.7)), 2)
})

test_that("output return type is correct",  {
  expect_equal(class(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)), "list")
  expect_equal(class(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$cov_c)[1], "matrix")
  expect_equal(class(gp_prior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)$cov_t)[1], "matrix")
  expect_equal(class(gp_posterior(Xc = xc, Xt = xt, Yc = yc, Yt = yt, sigma_hat = 1, choice = 1, l = 1)), "list")
  expect_equal(class(create_plot(X = x, Y = y, b = 1, col_num = 1, sigma_gp = 2, sigma_hat = 1.2, choice = 1, l = 0.7)), "list")
  expect_equal(class(create_plot(X = x, Y = y, b = 1, col_num = 1, sigma_gp = 2, sigma_hat = 1.2, choice = 1, l = 0.7))[1], "list")

})


test_that("expected errors", {
  expect_error(create_plot(X = x, Y = y, b = 1, col_num = 1, sigma_gp = 2, sigma_hat = 1.2, choice = 2, l = 1, alpha = 0), "Alpha cannot be zero")
  expect_error(create_plot(X = x, Y = y, b = 1, col_num = 1, sigma_gp = 2, sigma_hat = 1.2, choice = 3, l = 1, alpha = 1), "Choice should be either 1 or 2")
  expect_error(create_plot(X = x, Y = y, b = 1, col_num = 1, sigma_gp = 2, sigma_hat = 1.2, choice = 1, l = 0, alpha = 0.4), "Lengthscale should be positive")
  expect_error(create_plot(X = x, Y = y, b = 1, col_num = 1, sigma_gp = 2, sigma_hat = 1.2, choice = 2, l = 0, alpha = 1), "Lengthscale should be positive")
})

