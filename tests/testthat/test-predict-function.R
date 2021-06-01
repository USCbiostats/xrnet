library(bigmemory)

context("test predict function works correctly")

test_that("predict returns estimates for penalties already fit by xrnet object", {
  main_penalty <- define_penalty(0, user_penalty = c(2, 1, 0.05))
  external_penalty <- define_penalty(1, user_penalty = c(0.2, 0.1, 0.05))

  test_control <- xrnet_control(tolerance = 1e-15)

  xrnet_object <- xrnet(
    x = xtest,
    y = ytest,
    external = ztest,
    family = "gaussian",
    penalty_main = main_penalty,
    penalty_external = external_penalty,
    control = test_control
  )

  test_pred <- predict(xrnet_object, p = 1, pext = 0.05, type = "coefficients")
  test_coef <- coef(xrnet_object, p = 1, pext = 0.05)

  expect_identical(drop(test_pred$betas), xrnet_object$betas[, 2, 3])
  expect_identical(drop(test_coef$betas), xrnet_object$betas[, 2, 3])
  expect_identical(drop(test_pred$beta0), xrnet_object$beta0[2, 3])
  expect_identical(drop(test_pred$alphas), xrnet_object$alphas[, 2, 3])
  expect_identical(drop(test_pred$alpha0), xrnet_object$alpha0[2, 3])
})

test_that("predict returns right predictions for penalties already fit by xrnet object", {
  main_penalty <- define_penalty(0, user_penalty = c(2, 1, 0.05))
  external_penalty <- define_penalty(1, user_penalty = c(0.2, 0.1, 0.05))

  test_control <- xrnet_control(tolerance = 1e-15)

  xrnet_object <- xrnet(
    x = xtest,
    y = ytest,
    external = ztest,
    family = "gaussian",
    penalty_main = main_penalty,
    penalty_external = external_penalty,
    control = test_control
  )

  xtest_big <- as.big.matrix(xtest)

  predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[2, 2], xrnet_object$betas[, 2, 2])
  pred_xrnet <- predict(xrnet_object, p = 1, pext = 0.1, newdata = xtest)
  pred_xrnet_big <- predict(xrnet_object, p = 1, pext = 0.1, newdata = xtest_big)
  pred_xrnet_sparse <- predict(xrnet_object, p = 1, pext = 0.1, newdata = xsparse)
  expect_equivalent(pred_xrnet, predy)
  expect_equivalent(pred_xrnet_big, predy)
  expect_equivalent(pred_xrnet_sparse, predy)

  predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[2, 3], xrnet_object$betas[, 2, 3])
  pred_xrnet <- predict(xrnet_object, p = 1, pext = 0.05, newdata = xtest)
  pred_xrnet_big <- predict(xrnet_object, p = 1, pext = 0.05, newdata = xtest_big)
  pred_xrnet_sparse <- predict(xrnet_object, p = 1, pext = 0.05, newdata = xsparse)
  expect_equivalent(pred_xrnet, predy)
  expect_equivalent(pred_xrnet_big, predy)
  expect_equivalent(pred_xrnet_sparse, predy)

  predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[1, 3], xrnet_object$betas[, 1, 3])
  pred_xrnet <- predict(xrnet_object, p = 2, pext = 0.05, newdata = xtest)
  pred_xrnet_big <- predict(xrnet_object, p = 2, pext = 0.05, newdata = xtest_big)
  pred_xrnet_sparse <- predict(xrnet_object, p = 2, pext = 0.05, newdata = xsparse)
  expect_equivalent(pred_xrnet, predy)
  expect_equivalent(pred_xrnet_big, predy)
  expect_equivalent(pred_xrnet_sparse, predy)
})

test_that("predict returns right predictions for penalties already fit by xrnet object, no external data", {
  main_penalty <- define_penalty(penalty_type = 0, user_penalty = c(2, 1, 0.05))

  test_control <- xrnet_control(tolerance = 1e-15)

  xrnet_object <- xrnet(
    x = xtest,
    y = ytest,
    family = "gaussian",
    intercept = c(T, F),
    penalty_main = main_penalty,
    control = test_control
  )

  predy1 <- cbind(1, xtest) %*% c(xrnet_object$beta0[1, 1], xrnet_object$betas[, 1, 1])
  predy2 <- cbind(1, xtest) %*% c(xrnet_object$beta0[2, 1], xrnet_object$betas[, 2, 1])
  predy3 <- cbind(1, xtest) %*% c(xrnet_object$beta0[3, 1], xrnet_object$betas[, 3, 1])
  pred_xrnet <- predict(xrnet_object, p = c(0.05, 1, 2), newdata = xtest)
  expect_equivalent(pred_xrnet, cbind(predy1, predy2, predy3))
})

test_that("predict returns right estimates for penalties already fit by tune_xrnet object", {
  main_penalty <- define_penalty(0, user_penalty = c(2, 1, 0.05))
  external_penalty <- define_penalty(1, user_penalty = c(0.2, 0.1, 0.05))

  test_control <- xrnet_control(tolerance = 1e-15)

  xrnet_object <- tune_xrnet(
    x = xtest,
    y = ytest,
    external = ztest,
    family = "gaussian",
    penalty_main = main_penalty,
    penalty_external = external_penalty,
    control = test_control
  )

  test_pred <- predict(xrnet_object, p = "opt", pext = "opt", type = "coefficients")
  test_coef <- coef(xrnet_object, p = "opt", pext = "opt")

  optl1 <- which(xrnet_object$fitted_model$penalty == xrnet_object$opt_penalty)
  optl2 <- which(xrnet_object$fitted_model$penalty_ext == xrnet_object$opt_penalty_ext)

  expect_identical(drop(test_pred$betas), xrnet_object$fitted_model$betas[, optl1, optl2])
  expect_identical(drop(test_coef$betas), xrnet_object$fitted_model$betas[, optl1, optl2])
  expect_identical(drop(test_pred$beta0), xrnet_object$fitted_model$beta0[optl1, optl2])
  expect_identical(drop(test_pred$alphas), xrnet_object$fitted_model$alphas[, optl1, optl2])
  expect_identical(drop(test_pred$alpha0), xrnet_object$fitted_model$alpha0[optl1, optl2])
})
