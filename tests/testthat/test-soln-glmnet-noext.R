context("compare coefficients to glmnet when no external data (gaussian)")

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

# Ridge Regression #

test_that("x standardized, intercept", {

  # coef_glmnet <- glmnet(x = x, y = y, family = "gaussian", alpha = 0, thresh = 1e-15)

  penalty <- define_penalty(penalty_type = 0, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 1],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    b0_glmnet[1],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      penalty_main = penalty,
      control = test_control
    )$beta0[10, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, intercept", {

  # coef_glmnet <- glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, alpha = 0, thresh = 1e-15)

  penalty <- define_penalty(penalty_type = 0, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 2],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    b0_glmnet[2],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$beta0[10, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x standardized, NO intercept", {

  # coef_glmnet <- glmnet(x = x, y = y, family = "gaussian", intercept = FALSE, alpha = 0, thresh = 1e-15)

  penalty <- define_penalty(penalty_type = 0, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 3],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(F, F),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, NO intercept", {

  # coef_glmnet <- glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, intercept = FALSE, alpha = 0, thresh = 1e-15)

  penalty <- define_penalty(penalty_type = 0, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 4],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(F, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    b0_glmnet[4],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(F, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$beta0[10, 1] * sd_y,
    tolerance = 1e-5
  )
})

# Lasso Regression #

test_that("x standardized, intercept", {

  # coef_glmnet <- glmnet(x = x, y = y, family = "gaussian", alpha = 1, thresh = 1e-15)

  penalty <- define_penalty(penalty_type = 1, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 5],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    b0_glmnet[5],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      penalty_main = penalty,
      control = test_control
    )$beta0[10, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, intercept", {

  # coef_glmnet <- glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, alpha = 1, thresh = 1e-15)

  penalty <- define_penalty(penalty_type = 1, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 6],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    b0_glmnet[6],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$beta0[10, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x standardized, NO intercept", {

  # coef_glmnet <- glmnet(x = x, y = y, family = "gaussian", intercept = FALSE, alpha = 1, thresh = 1e-15)

  penalty <- define_penalty(penalty_type = 1, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 7],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(F, F),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, NO intercept", {

  # coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, intercept = FALSE, alpha = 1, thresh = 1e-15)

  penalty <- define_penalty(penalty_type = 1, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 8],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(F, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )
})

# Elastic Net Regression #

test_that("x standardized, intercept", {

  # coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", alpha = 0.5, thresh = 1e-15)
  # betas_glmnet9 <- drop(unname(as.matrix(coef_glmnet$beta[,10])))
  # b0_glmnet9 <- unname(coef_glmnet$a0[10])

  penalty <- define_penalty(penalty_type = 0.5, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 9],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    b0_glmnet[9],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      penalty_main = penalty,
      control = test_control
    )$beta0[10, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, intercept", {

  # coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, alpha = 0.5, thresh = 1e-15)
  # betas_glmnet10 <- drop(unname(as.matrix(coef_glmnet$beta[,10])))
  # b0_glmnet10 <- unname(coef_glmnet$a0[10])

  penalty <- define_penalty(penalty_type = 0.5, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 10],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )

  expect_equal(
    b0_glmnet[10],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(T, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$beta0[10, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x standardized, NO intercept", {

  # coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", intercept = FALSE, alpha = 0.5, thresh = 1e-15)
  # betas_glmnet11 <- drop(unname(as.matrix(coef_glmnet$beta[,10])))
  # b0_glmnet11 <- unname(coef_glmnet$a0[10])

  penalty <- define_penalty(penalty_type = 0.5, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 11],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(F, F),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, NO intercept", {

  # coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, intercept = FALSE, alpha = 0.5, thresh = 1e-15)
  # betas_glmnet12 <- drop(unname(as.matrix(coef_glmnet$beta[,10])))
  # b0_glmnet12 <- unname(coef_glmnet$a0[10])

  penalty <- define_penalty(penalty_type = 0.5, num_penalty = 100)

  expect_equal(
    betas_glmnet[, 12],
    xrnet(
      x = xtest,
      y = ytest_scaled,
      family = "gaussian",
      intercept = c(F, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1] * sd_y,
    tolerance = 1e-5
  )
})


# Ridge - No Penalty on 1st two variables #

test_that("x NOT standardized, intercept", {
  pf <- c(rep(0, 2), rep(1, NCOL(xtest) - 2))
  # coef_glmnet <- glmnet(x = xtest, y = ytest, family = "gaussian", standardize = FALSE, alpha = 0, thresh = 1e-15, penalty.factor = pf)

  penalty <- define_penalty(
    penalty_type = 0,
    user_penalty = path_ridge,
    custom_multiplier = NCOL(xtest) * pf[-c(1, 2)] / sum(pf)
  )

  expect_equal(
    betas_glmnet[-c(1, 2), 13],
    xrnet(
      x = xtest[, -c(1, 2)],
      y = ytest,
      unpen = xtest[, c(1, 2)],
      family = "gaussian",
      intercept = c(T, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1],
    tolerance = 1e-5
  )
})

# Lasso - No Penalty on 1st two variables #

test_that("x NOT standardized, intercept", {
  pf <- c(rep(0, 2), rep(1, NCOL(xtest) - 2))
  # coef_glmnet <- glmnet(x = xtest, y = ytest, family = "gaussian", standardize = FALSE, alpha = 1, thresh = 1e-15,penalty.factor = pf)

  penalty <- define_penalty(
    penalty_type = 1,
    user_penalty = path_lasso,
    custom_multiplier = NCOL(xtest) * pf[-c(1, 2)] / sum(pf)
  )

  expect_equal(
    betas_glmnet[-c(1, 2), 14],
    xrnet(
      x = xtest[, -c(1, 2)],
      y = ytest,
      unpen = xtest[, c(1, 2)],
      family = "gaussian",
      intercept = c(T, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1],
    tolerance = 1e-5
  )
})

# Elastic Net - No Penalty on 1st two variables #

test_that("x NOT standardized, intercept", {
  pf <- c(rep(0, 2), rep(1, NCOL(xtest) - 2))
  # coef_glmnet <- glmnet(x = xtest, y = ytest, family = "gaussian", standardize = FALSE, alpha = 0.5, thresh = 1e-15, penalty.factor = pf)

  penalty <- define_penalty(
    penalty_type = 0.5,
    user_penalty = path_en,
    custom_multiplier = NCOL(xtest) * pf[-c(1, 2)] / sum(pf)
  )

  expect_equal(
    betas_glmnet[-c(1, 2), 15],
    xrnet(
      x = xtest[, -c(1, 2)],
      y = ytest,
      unpen = xtest[, c(1, 2)],
      family = "gaussian",
      intercept = c(T, F),
      standardize = c(F, T),
      penalty_main = penalty,
      control = test_control
    )$betas[, 10, 1],
    tolerance = 1e-5
  )
})
