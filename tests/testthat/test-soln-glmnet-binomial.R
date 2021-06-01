library(Matrix)

# Code to generate glmnet solutions
#
# library(glmnet)
#
# fit_glmnet1 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 0,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet2 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 0,
#     standardize = FALSE,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet3 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 0,
#     intercept = FALSE,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet4 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 0,
#     intercept = FALSE,
#     standardize = FALSE,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet5 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 1,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet6 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 1,
#     standardize = FALSE,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet7 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 1,
#     intercept = FALSE,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet8 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 1,
#     intercept = FALSE,
#     standardize = FALSE,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet9 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 0.5,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet10 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 0.5,
#     standardize = FALSE,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet11 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 0.5,
#     intercept = FALSE,
#     lambda.min.ratio = 0.01
# )
#
# fit_glmnet12 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 0.5,
#     intercept = FALSE,
#     standardize = FALSE,
#     lambda.min.ratio = 0.01
# )
#
# pf <- c(rep(0, 2), rep(1, NCOL(xtest_binomial) - 2))
#
# fit_glmnet13 <- glmnet(
#     x = xtest_binomial,
#     y = ytest_binomial,
#     family = "binomial",
#     thresh = 1e-15,
#     alpha = 0.5,
#     penalty.factor = pf,
#     lambda.min.ratio = 0.01
# )
#
# betas_binomial <- matrix(NA, nrow = NCOL(xtest_binomial) + 1, ncol = 13)
#
# for (i in 1:13) {
#     fit_current <- get(paste0("fit_glmnet", i))
#     betas_binomial[, i] <- as.vector(coef(fit_current, s = fit_current$lambda[10]))
# }
#
# saveRDS(betas_binomial, "tests/testthat/testdata/betas_binomial.rds")
# saveRDS(fit_glmnet13$lambda, "tests/testthat/testdata/lambda_binomial.rds")


context("compare coefficients to glmnet when no external data (binomial)")

# Ridge Regression #

test_that("x standardized, intercept, ridge", {
  penalty <- define_ridge(num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    control = xrnet_control(tolerance = 1e-15)
  )

  fit_xrnet_sparse <- xrnet(
    x = Matrix(xtest_binomial, sparse = TRUE),
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 1],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 1],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[1, 1],
    drop(fit_xrnet_sparse$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 1],
    drop(fit_xrnet_sparse$betas)[, 10],
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, intercept, ridge", {
  penalty <- define_ridge(num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    standardize = c(FALSE, FALSE),
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 2],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 2],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

test_that("x standardized, NO intercept, ridge", {
  penalty <- define_ridge(num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    intercept = c(FALSE, FALSE),
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 3],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 3],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, NO intercept, ridge", {
  penalty <- define_ridge(num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    standardize = c(FALSE, FALSE),
    intercept = c(FALSE, FALSE),
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 4],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 4],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

# Lasso Regression #

test_that("x standardized, intercept, lasso", {
  penalty <- define_lasso(num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 5],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 5],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, intercept, lasso", {
  penalty <- define_lasso(num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    standardize = c(FALSE, FALSE),
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 6],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 6],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

test_that("x standardized, NO intercept, lasso", {
  penalty <- define_lasso(num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    intercept = c(FALSE, FALSE),
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 7],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 7],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, NO intercept, lasso", {
  penalty <- define_lasso(num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    standardize = c(FALSE, FALSE),
    intercept = c(FALSE, FALSE),
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 8],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 8],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

# Elastic Net Regression #

test_that("x standardized, intercept, en", {
  penalty <- define_enet(en_param = 0.5, num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 9],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 9],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, intercept, en", {
  penalty <- define_enet(en_param = 0.5, num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    standardize = c(FALSE, FALSE),
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 10],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 10],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

test_that("x standardized, NO intercept, en", {
  penalty <- define_enet(en_param = 0.5, num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    intercept = c(FALSE, FALSE),
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 11],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 11],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

test_that("x NOT standardized, NO intercept, en", {
  penalty <- define_penalty(penalty_type = 0.5, num_penalty = 100, penalty_ratio = 0.01)

  fit_xrnet <- xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = penalty,
    standardize = c(FALSE, FALSE),
    intercept = c(FALSE, FALSE),
    control = xrnet_control(tolerance = 1e-15)
  )

  expect_equal(
    betas_binomial[1, 12],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 12],
    drop(fit_xrnet$betas)[, 10],
    tolerance = 1e-5
  )
})

# Elastic Net - No Penalty on 1st two variables #

test_that("x standardized, intercept, unpenalized variables", {
  pf <- c(rep(0, 2), rep(1, NCOL(xtest_binomial) - 2))

  penalty <- define_enet(
    en_param = 0.5,
    user_penalty = lambda_binomial,
    custom_multiplier = (NCOL(xtest_binomial) / sum(pf)) * rep(1, NCOL(xtest_binomial) - 2)
  )

  fit_xrnet <- xrnet(
    x = xtest_binomial[, -c(1, 2)],
    y = ytest_binomial,
    unpen = xtest_binomial[, c(1, 2)],
    family = "binomial",
    penalty_main = penalty,
    control = xrnet_control(tolerance = 1e-15)
  )

  betas_all <- c(drop(fit_xrnet$gammas)[, 10], drop(fit_xrnet$betas)[, 10])

  expect_equal(
    betas_binomial[1, 13],
    drop(fit_xrnet$beta0)[10],
    tolerance = 1e-5
  )

  expect_equal(
    betas_binomial[-1, 13],
    betas_all,
    tolerance = 1e-5
  )
})
