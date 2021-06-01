context("compare coefficent estimates to CVX (manual penalty)")

##### Code used to generate data files for all tests #####

# set.seed(123)
# n <- 100
# p <- 50
# q <- 5

# meanx_true <- rep(0, p)
# covx_true <- matrix(NA, nrow = p, ncol = p)
# for (i in 1:p) {
#     for (j in 1:p) {
#         covx_true[i, j] <- 0.5^abs(i-j)
#     }
# }

# meanz_true <- rep(0, q)
# covz_true <- matrix(NA, nrow = q, ncol = q)
# for (i in 1:q) {
#     for (j in 1:q) {
#         covz_true[i, j] <- 0.5^abs(i-j)
#     }
# }

# a0 <- 0.01
# a <- c(0.1, -0.1, rep(0, q - 2))
# z <- mvrnorm(n = p, mu = meanz_true, Sigma = covz_true)
# b <- drop(a0 + z %*% a + 0.2*rnorm(p))
# x <- mvrnorm(n = n, mu = meanx_true, Sigma = covx_true)
# y <- drop(x %*% b + rnorm(n))

# mean_x <- colMeans(x)
# var_x <- apply(x, 2, var) * (n-1) / n
# sd_x <- sqrt(var_x)

# xscaled <- matrix(NA, nrow = n, ncol = p)
# for (i in 1:p) {
#     xscaled[, i] <- (x[, i] - mean_x[i]) / sd_x[i]
# }

# mean_z <- colMeans(z)
# sd_z <- sqrt(apply(z, 2, var) * (p-1) / p)

# zscaled <- matrix(NA, nrow = p, ncol = q)
# for (i in 1:q) {
#     zscaled[, i] <- (z[, i] - mean_z[i]) / sd_z[i]
# }

test_that("x and ext standardized, both intercepts", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(
    alphas_cvx_mat[, 1],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      intercept = c(T, T),
      standardize = c(T, T),
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 1],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      intercept = c(T, T),
      standardize = c(T, T),
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    alphas_cvx_mat[, 1],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      intercept = c(T, T),
      standardize = c(T, T),
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 1],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      intercept = c(T, T),
      standardize = c(T, T),
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, ext standardized, both intercepts", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(
    alphas_cvx_mat[, 2],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(F, T),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 2],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(F, T),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    alphas_cvx_mat[, 2],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(F, T),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 2],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(F, T),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x standardized, ext NOT standardized, both intercepts", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(
    alphas_cvx_mat[, 3],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(T, F),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 3],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(T, F),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    alphas_cvx_mat[, 3],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(T, F),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 3],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(T, F),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, ext NOT standardized, both intercepts", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(
    alphas_cvx_mat[, 4],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(F, F),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 4],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(F, F),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    alphas_cvx_mat[, 4],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(F, F),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 4],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(T, T),
      standardize = c(F, F),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x standardized, ext standardized, no 2nd level intercept", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(
    alphas_cvx_mat[, 5],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(TRUE, TRUE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 5],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(TRUE, TRUE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    alphas_cvx_mat[, 5],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(TRUE, TRUE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 5],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(TRUE, TRUE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, ext standardized, no 2nd level intercept", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(
    alphas_cvx_mat[, 6],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(FALSE, TRUE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 6],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(FALSE, TRUE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    alphas_cvx_mat[, 6],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(FALSE, TRUE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 6],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(FALSE, TRUE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x standardized, ext NOT standardized, no 2nd level intercept", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(
    alphas_cvx_mat[, 7],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(TRUE, FALSE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 7],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(TRUE, FALSE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    alphas_cvx_mat[, 7],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(TRUE, FALSE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 7],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(TRUE, FALSE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, ext NOT standardized, no 2nd level intercept", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(
    alphas_cvx_mat[, 8],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(FALSE, FALSE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 8],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(FALSE, FALSE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    alphas_cvx_mat[, 8],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(FALSE, FALSE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    betas_cvx_mat[, 8],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(TRUE, FALSE),
      standardize = c(FALSE, FALSE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x standardized, ext standardized, no 1st level intercept", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(alphas_cvx_mat[, 9],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(TRUE, TRUE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(betas_cvx_mat[, 9],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(TRUE, TRUE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(alphas_cvx_mat[, 9],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(TRUE, TRUE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(betas_cvx_mat[, 9],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(TRUE, TRUE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, ext standardized, no 1st level intercept", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(alphas_cvx_mat[, 10],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(FALSE, TRUE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(betas_cvx_mat[, 10],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(FALSE, TRUE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(alphas_cvx_mat[, 10],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(FALSE, TRUE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(betas_cvx_mat[, 10],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(FALSE, TRUE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x standardized, ext NOT standardized, no 1st level intercept", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(alphas_cvx_mat[, 11],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(TRUE, FALSE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(betas_cvx_mat[, 11],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(TRUE, FALSE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(alphas_cvx_mat[, 11],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(TRUE, FALSE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(betas_cvx_mat[, 11],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(TRUE, FALSE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, ext NOT standardized, no 1st level intercept", {
  test_control <- list(tolerance = 1e-20)

  expect_equal(alphas_cvx_mat[, 12],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(FALSE, FALSE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(betas_cvx_mat[, 12],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(FALSE, FALSE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(alphas_cvx_mat[, 12],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(FALSE, FALSE),
      control = test_control
    )$alphas[1:5, 1, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(betas_cvx_mat[, 12],
    xrnet(
      x = xsparse,
      y = ytest_scaled,
      external = ztest,
      family = "gaussian",
      penalty_main = define_penalty(0, user_penalty = 1),
      penalty_external = define_penalty(1, user_penalty = 0.1),
      intercept = c(FALSE, TRUE),
      standardize = c(FALSE, FALSE),
      control = test_control
    )$betas[1:50, 1, 1] * sd_y,
    tolerance = 1e-5
  )
})
