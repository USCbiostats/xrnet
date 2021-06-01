library(doParallel)
library(bigmemory)

context("check computation of CV fold errors")

test_that("gaussian, mse (sequential)", {
  main_penalty <- define_penalty(0, num_penalty = 20)
  external_penalty <- define_penalty(1, num_penalty = 20)

  fit_xrnet <- tune_xrnet(
    x = xtest,
    y = ytest,
    external = ztest,
    family = "gaussian",
    penalty_main = main_penalty,
    penalty_external = external_penalty,
    control = list(tolerance = 1e-10),
    loss = "mse",
    foldid = foldid
  )

  expect_equal(
    which.min(cv_mean),
    which.min(fit_xrnet$cv_mean),
    check.attribute = FALSE
  )

  fit_xrnet <- tune_xrnet(
    x = xsparse,
    y = ytest,
    external = ztest,
    family = "gaussian",
    penalty_main = main_penalty,
    penalty_external = external_penalty,
    control = list(tolerance = 1e-10),
    loss = "mse",
    foldid = foldid
  )

  expect_equal(
    which.min(cv_mean),
    which.min(fit_xrnet$cv_mean),
    check.attribute = FALSE
  )

  fit_xrnet <- tune_xrnet(
    x = xtest,
    y = ytest,
    external = zsparse,
    family = "gaussian",
    penalty_main = main_penalty,
    penalty_external = external_penalty,
    control = list(tolerance = 1e-10),
    loss = "mse",
    foldid = foldid
  )

  expect_equal(
    which.min(cv_mean),
    which.min(fit_xrnet$cv_mean),
    check.attribute = FALSE
  )
})

test_that("gaussian, mse (parallel)", {
  main_penalty <- define_penalty(0, num_penalty = 20)
  external_penalty <- define_penalty(1, num_penalty = 20)

  cl <- makeCluster(2, type = "PSOCK")
  registerDoParallel(cl)

  fit_xrnet <- tune_xrnet(
    x = xtest,
    y = ytest,
    external = ztest,
    family = "gaussian",
    penalty_main = main_penalty,
    penalty_external = external_penalty,
    control = list(tolerance = 1e-10),
    loss = "mse",
    foldid = foldid,
    parallel = TRUE
  )

  expect_equal(
    which.min(cv_mean),
    which.min(fit_xrnet$cv_mean),
    check.attribute = FALSE
  )

  fit_xrnet <- tune_xrnet(
    x = as.big.matrix(xtest),
    y = ytest,
    external = ztest,
    family = "gaussian",
    penalty_main = main_penalty,
    penalty_external = external_penalty,
    control = list(tolerance = 1e-10),
    loss = "mse",
    foldid = foldid,
    parallel = TRUE
  )

  expect_equal(
    which.min(cv_mean),
    which.min(fit_xrnet$cv_mean),
    check.attribute = FALSE
  )
})

test_that("gaussian, mae (sequential)", {
  main_penalty <- define_penalty(0, num_penalty = 20)

  fit_xrnet <- tune_xrnet(
    x = xtest,
    y = ytest,
    family = "gaussian",
    penalty_main = main_penalty,
    control = list(tolerance = 1e-12),
    loss = "mae",
    foldid = foldid
  )

  expect_equal(
    cv_mae,
    drop(fit_xrnet$cv_mean),
    check.attribute = FALSE,
    tolerance = 1e-5
  )
})

test_that("binomial, auc (sequential)", {
  main_penalty <- define_penalty(
    penalty_type = 0,
    num_penalty = 20,
    penalty_ratio = 0.001
  )

  fit_xrnet <- tune_xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = main_penalty,
    control = list(tolerance = 1e-10),
    loss = "auc",
    foldid = foldid_binomial
  )

  expect_equal(
    cv_auc,
    drop(fit_xrnet$cv_mean),
    check.attribute = FALSE,
    tolerance = 1e-5
  )
})

test_that("binomial, deviance (sequential)", {
  main_penalty <- define_penalty(
    penalty_type = 0,
    num_penalty = 20,
    penalty_ratio = 0.001
  )

  fit_xrnet <- tune_xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = main_penalty,
    control = list(tolerance = 1e-10),
    loss = "deviance",
    foldid = foldid_binomial
  )

  expect_equal(
    cv_deviance,
    drop(fit_xrnet$cv_mean),
    check.attribute = FALSE
  )
})
