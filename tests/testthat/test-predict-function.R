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
    skip('internal')
    myPenalty <- define_penalty(penalty_type = 0,
                                penalty_type_ext = 1,
                                user_penalty = c(1:20),
                                user_penalty_ext = seq(0.01, 0.1, 0.01))

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

    test <- predict(xrnet_object2,
                    p = 4,
                    pext = 0.05,
                    type = "coefficients",
                    penalty = myPenalty2,
                    x = xtest,
                    y = ytest,
                    exteranl = ztest,
                    control = myControl)

    expect_equal(drop(predict(xrnet_object2,
                              p = 4,
                              pext = 0.05,
                              type = "coefficients",
                              penalty = myPenalty,
                              x = xtest,
                              y = ytest,
                              external = ztest,
                              family = "gaussian",
                              control = myControl)$betas),
                 xrnet_object$betas[, 17, 6],
                 tolerance = 1e-6)

    expect_equal(drop(coef(xrnet_object2,
                           p = 4,
                           pext = 0.05,
                           penalty = myPenalty,
                           x = xtest,
                           y = ytest,
                           external = ztest,
                           family = "gaussian",
                           control = myControl)$betas),
                 xrnet_object$betas[, 17, 6],
                 tolerance = 1e-6)

    expect_equal(drop(test$beta0), xrnet_object$beta0[17, 6], tolerance = 1e-6)
    expect_equal(drop(test$alphas), xrnet_object$alphas[ , 17, 6], tolerance = 1e-6)
    expect_equal(drop(test$alpha0), xrnet_object$alpha0[17, 6], tolerance = 1e-6)
})

test_that("predict returns right predictions for penalties not already fit by xrnet object", {
    skip('internal')
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

    myPenalty2 <- define_penalty(penalty_type = 0,
                                 penalty_type_ext = 1,
                                 user_penalty = c(4, 3),
                                 user_penalty_ext = c(0.25, 0.15))

    xrnet_object2 <- xrnet(x = xtest,
                           y = ytest,
                           external = ztest,
                           family = "gaussian",
                           penalty = myPenalty2,
                           control = myControl)

    predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[2, 2], xrnet_object$betas[, 2, 2])
    expect_equivalent(predict(xrnet_object2,
                              p = 1,
                              pext = 0.1,
                              newdata = xtest,
                              x = xtest,
                              y = ytest,
                              external = ztest,
                              penalty = myPenalty,
                              control = myControl),
                      predy,
                      tolerance = 1e-6)

    predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[2, 3], xrnet_object$betas[, 2, 3])

    expect_equivalent(predict(xrnet_object2,
                              p = 1,
                              pext = 0.05,
                              newdata = xtest,
                              x = xtest,
                              y = ytest,
                              external = ztest,
                              penalty = myPenalty,
                              control = myControl),
                      predy,
                      tolerance = 1e-6)

    predy <- cbind(1, xtest) %*% c(xrnet_object$beta0[1, 3], xrnet_object$betas[, 1, 3])
    expect_equivalent(predict(xrnet_object2,
                              p = 2,
                              pext = 0.05,
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
