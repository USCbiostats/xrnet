context("test predict function works correctly")

test_that("predict returns estimates for penalties already fit by xrnet object", {
    myPenalty <- define_penalty(penalty_type = 0,
                                penalty_type_ext = 1,
                                user_penalty = c(2, 1, 0.05),
                                user_penalty_ext = c(0.2, 0.1, 0.05))

    myControl <- xrnet.control(tolerance = 1e-15)

    xrnet_object <- xrnet(x = xtest,
                          y = ytest,
                          external = ztest,
                          family = "gaussian",
                          penalty = myPenalty,
                          control = myControl)

    test <- predict(xrnet_object, p = 1, pext = 0.05, type = "coefficients")
    expect_identical(drop(predict(xrnet_object, p = 1, pext = 0.05, type = "coefficients")$betas), xrnet_object$betas[, 2, 3])
    expect_identical(drop(coef(xrnet_object, p = 1, pext = 0.05)$betas), xrnet_object$betas[, 2, 3])
    expect_identical(drop(test$beta0), xrnet_object$beta0[2, 3])
    expect_identical(drop(test$alphas), xrnet_object$alphas[ , 2, 3])
    expect_identical(drop(test$alpha0), xrnet_object$alpha0[2, 3])
})

test_that("predict returns right predictions for penalties already fit by xrnet object", {
    myPenalty <- define_penalty(penalty_type = 0,
                                penalty_type_ext = 1,
                                user_penalty = c(2, 1, 0.05),
                                user_penalty_ext = c(0.2, 0.1, 0.05))

    myControl <- xrnet.control(tolerance = 1e-15)

    xrnet_object <- xrnet(x = xtest,
                          y = ytest,
                          external = ztest,
                          family = "gaussian",
                          penalty = myPenalty,
                          control = myControl)

    predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[2, 2], xrnet_object$betas[, 2, 2])
    expect_equivalent(predict(xrnet_object, p = 1, pext = 0.1, newdata = xtest), predy)

    predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[2, 3], xrnet_object$betas[, 2, 3])
    expect_equivalent(predict(xrnet_object, p = 1, pext = 0.05, newdata = xtest), predy)

    predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[1, 3], xrnet_object$betas[, 1, 3])
    expect_equivalent(predict(xrnet_object, p = 2, pext = 0.05, newdata = xtest), predy)
})

test_that("predict returns estimates for penalties not fit by xrnet object", {
    skip("Internal")
    myPenalty <- define_penalty(penalty_type = 0,
                                penalty_type_ext = 1,
                                num_penalty = 25,
                                num_penalty_ext = 25)

    myControl <- xrnet.control(tolerance = 1e-15)

    xrnet_object <- xrnet(x = xtest,
                          y = ytest,
                          external = ztest,
                          family = "gaussian",
                          penalty = myPenalty,
                          control = myControl)

    myPenalty2 <- define_penalty()

    xrnet_object2 <- xrnet(x = xtest,
                           y = ytest,
                           external = ztest,
                           family = "gaussian",
                           penalty = myPenalty2,
                           control = myControl)

    idx1 <- 10
    idx2 <- 15

    test <- predict(xrnet_object2,
                    p = xrnet_object$penalty[idx1],
                    pext = xrnet_object$penalty_ext[idx2],
                    type = "coefficients",
                    penalty = myPenalty2,
                    x = xtest,
                    y = ytest,
                    external = ztest,
                    control = myControl)

    expect_equal(drop(predict(xrnet_object2,
                              p = xrnet_object$penalty[idx1],
                              pext = xrnet_object$penalty_ext[idx2],
                              type = "coefficients",
                              penalty = myPenalty,
                              x = xtest,
                              y = ytest,
                              external = ztest,
                              family = "gaussian",
                              control = myControl)$betas),
                 xrnet_object$betas[, 10, 15],
                 tolerance = 1e-6)

    expect_equal(drop(coef(xrnet_object2,
                           p = xrnet_object$penalty[idx1],
                           pext = xrnet_object$penalty_ext[idx2],
                           penalty = myPenalty,
                           x = xtest,
                           y = ytest,
                           external = ztest,
                           family = "gaussian",
                           control = myControl)$betas),
                 xrnet_object$betas[, idx1, idx2],
                 tolerance = 1e-6)

    expect_equal(drop(test$beta0), xrnet_object$beta0[idx1, idx2], tolerance = 1e-6)
    expect_equal(drop(test$alphas), xrnet_object$alphas[ , idx1, idx2], tolerance = 1e-6)
    expect_equal(drop(test$alpha0), xrnet_object$alpha0[idx1, idx2], tolerance = 1e-6)
})

test_that("predict returns right predictions for penalties not already fit by xrnet object", {
    skip("Internal")
    myPenalty <- define_penalty(penalty_type = 0,
                                penalty_type_ext = 1,
                                num_penalty = 15,
                                num_penalty_ext = 15)

    myControl <- xrnet.control(tolerance = 1e-20)

    xrnet_object <- xrnet(x = xtest,
                          y = ytest,
                          external = ztest,
                          family = "gaussian",
                          penalty = myPenalty,
                          control = myControl)

    myPenalty2 <- define_penalty(penalty_type = 0,
                                 penalty_type_ext = 1,
                                 num_penalty = 20,
                                 num_penalty_ext = 20)

    xrnet_object2 <- xrnet(x = xtest,
                           y = ytest,
                           external = ztest,
                           family = "gaussian",
                           penalty = myPenalty2,
                           control = myControl)

    idx1 <- 3
    idx2 <- 2

    predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[idx1, idx2], xrnet_object$betas[, idx1, idx2])
    expect_equivalent(predict(xrnet_object2,
                              p = xrnet_object$penalty[idx1],
                              pext = xrnet_object$penalty_ext[idx2],
                              newdata = xtest,
                              x = xtest,
                              y = ytest,
                              external = ztest,
                              penalty = myPenalty,
                              control = myControl),
                      predy,
                      tolerance = 1e-6)

    predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[idx1 + 1, idx2 + 1], xrnet_object$betas[, idx1 + 1, idx2 + 1])
    expect_equivalent(predict(xrnet_object2,
                              p = xrnet_object$penalty[idx1 + 1],
                              pext = xrnet_object$penalty_ext[idx2 + 1],
                              newdata = xtest,
                              x = xtest,
                              y = ytest,
                              external = ztest,
                              penalty = myPenalty,
                              control = myControl),
                      predy,
                      tolerance = 1e-6)

    predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[idx1 + 2, idx2 + 2], xrnet_object$betas[, idx1 + 2, idx2 + 2])
    expect_equivalent(predict(xrnet_object2,
                              p = xrnet_object$penalty[idx1 + 2],
                              pext = xrnet_object$penalty_ext[idx2 + 2],
                              newdata = xtest,
                              x = xtest,
                              y = ytest,
                              external = ztest,
                              penalty = myPenalty,
                              control = myControl),
                      predy,
                      tolerance = 1e-6)
})

test_that("predict returns right predictions for penalties already fit by xrnet object, no external data", {
    myPenalty <- define_penalty(penalty_type = 0,
                                user_penalty = c(2, 1, 0.05))

    myControl <- xrnet.control(tolerance = 1e-15)

    xrnet_object <- xrnet(x = xtest,
                          y = ytest,
                          family = "gaussian",
                          intercept = c(T, F),
                          penalty = myPenalty,
                          control = myControl)

    predy1 <- cbind(1, xtest) %*% c(xrnet_object$beta0[1, 1], xrnet_object$betas[, 1, 1])
    predy2 <- cbind(1, xtest) %*% c(xrnet_object$beta0[2, 1], xrnet_object$betas[, 2, 1])
    predy3 <- cbind(1, xtest) %*% c(xrnet_object$beta0[3, 1], xrnet_object$betas[, 3, 1])
    expect_equivalent(
        predict(xrnet_object,
                p = c(0.05, 1, 2),
                newdata = xtest),
        cbind(predy1, predy2, predy3)
    )
})

test_that("predict returns estimates for penalties already fit by tune_xrnet object", {
    myPenalty <- define_penalty(
        penalty_type = 0,
        penalty_type_ext = 1,
        user_penalty = c(2, 1, 0.05),
        user_penalty_ext = c(0.2, 0.1, 0.05)
    )

    myControl <- xrnet.control(tolerance = 1e-15)

    xrnet_object <- tune_xrnet(
        x = xtest,
        y = ytest,
        external = ztest,
        family = "gaussian",
        penalty = myPenalty,
        control = myControl
    )

    test <- predict(xrnet_object, p = 1, pext = 0.05, type = "coefficients")
    expect_identical(drop(predict(xrnet_object, p = 1, pext = 0.05, type = "coefficients")$betas), xrnet_object$fitted_model$betas[, 2, 3])
    expect_identical(drop(coef(xrnet_object, p = 1, pext = 0.05)$betas), xrnet_object$fitted_model$betas[, 2, 3])
    expect_identical(drop(test$beta0), xrnet_object$fitted_model$beta0[2, 3])
    expect_identical(drop(test$alphas), xrnet_object$fitted_model$alphas[ , 2, 3])
    expect_identical(drop(test$alpha0), xrnet_object$fitted_model$alpha0[2, 3])
})
