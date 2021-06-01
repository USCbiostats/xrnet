library(doParallel)
library(bigmemory)

context("check computation of CV fold errors")

test_that("gaussian, mse (sequential)", {
  mainPen <- define_penalty(0, num_penalty = 20)
  extPen <- define_penalty(1, num_penalty = 20)

  fit_xrnet <- tune_xrnet(
    x = xtest,
    y = ytest,
    external = ztest,
    family = "gaussian",
    penalty_main = mainPen,
    penalty_external = extPen,
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
    penalty_main = mainPen,
    penalty_external = extPen,
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
    penalty_main = mainPen,
    penalty_external = extPen,
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
  mainPen <- define_penalty(0, num_penalty = 20)
  extPen <- define_penalty(1, num_penalty = 20)

  cl <- makeCluster(2, type = "PSOCK")
  registerDoParallel(cl)

  fit_xrnet <- tune_xrnet(
    x = xtest,
    y = ytest,
    external = ztest,
    family = "gaussian",
    penalty_main = mainPen,
    penalty_external = extPen,
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
    penalty_main = mainPen,
    penalty_external = extPen,
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
  mainPen <- define_penalty(0, num_penalty = 20)

  fit_xrnet <- tune_xrnet(
    x = xtest,
    y = ytest,
    family = "gaussian",
    penalty_main = mainPen,
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
  mainPen <- define_penalty(
    penalty_type = 0,
    num_penalty = 20,
    penalty_ratio = 0.001
  )

  fit_xrnet <- tune_xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = mainPen,
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
  mainPen <- define_penalty(
    penalty_type = 0,
    num_penalty = 20,
    penalty_ratio = 0.001
  )

  fit_xrnet <- tune_xrnet(
    x = xtest_binomial,
    y = ytest_binomial,
    family = "binomial",
    penalty_main = mainPen,
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
