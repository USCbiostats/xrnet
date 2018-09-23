context("test predict function works correctly")

xtest <- readRDS(file = "testdata/xtest.rds")
ytest <- readRDS(file = "testdata/ytest.rds")
ztest <- readRDS(file = "testdata/ztest.rds")

myControl <- list(tolerance = 1e-15, earlyStop = FALSE)

test_that("predict returns estimates for penalties already fit by hierr object", {
    myPenalty <- definePenalty(penalty_type = 0,
                               penalty_type_ext = 1,
                               user_penalty = c(2, 1, 0.05),
                               user_penalty_ext = c(0.2, 0.1, 0.05))
    hierr_object <- hierr(x = xtest,
                          y = ytest,
                          external = ztest,
                          family = "gaussian",
                          penalty = myPenalty,
                          control = myControl)

    test <- predict(hierr_object, p = 1, pext = 0.05, type = "coefficients")
    expect_identical(drop(predict(hierr_object, p = 1, pext = 0.05, type = "coefficients")$betas), hierr_object$betas[, 2, 3])
    expect_identical(drop(test$beta0), hierr_object$beta0[2, 3])
    expect_identical(drop(test$alphas), hierr_object$alphas[ , 2, 3])
    expect_identical(drop(test$alpha0), hierr_object$alpha0[2, 3])
})

test_that("predict returns right predictions for penalties already fit by hierr object", {
    myPenalty <- definePenalty(penalty_type = 0,
                               penalty_type_ext = 1,
                               user_penalty = c(2, 1, 0.05),
                               user_penalty_ext = c(0.2, 0.1, 0.05))

    hierr_object <- hierr(x = xtest,
                          y = ytest,
                          external = ztest,
                          family = "gaussian",
                          penalty = myPenalty,
                          control = myControl)

    predy <- cbind(1, xtest) %*% c(hierr_object$beta0[2, 2], hierr_object$betas[, 2, 2])
    expect_equivalent(predict(hierr_object, p = 1, pext = 0.1, newdata = xtest), predy)

    predy <- cbind(1, xtest) %*% c(hierr_object$beta0[2, 3], hierr_object$betas[, 2, 3])
    expect_equivalent(predict(hierr_object, p = 1, pext = 0.05, newdata = xtest), predy)

    predy <- cbind(1, xtest) %*% c(hierr_object$beta0[1, 3], hierr_object$betas[, 1, 3])
    expect_equivalent(predict(hierr_object, p = 2, pext = 0.05, newdata = xtest), predy)
})

test_that("predict returns estimates for penalties not fit by hierr object", {
    myPenalty <- definePenalty(penalty_type = 0,
                               penalty_type_ext = 1,
                               user_penalty = c(1:20),
                               user_penalty_ext = seq(0.01, 0.1, 0.01))

    hierr_object <- hierr(x = xtest,
                          y = ytest,
                          external = ztest,
                          family = "gaussian",
                          penalty = myPenalty,
                          control = myControl)

    hierr_object2 <- hierr(x = xtest,
                           y = ytest,
                           external = ztest,
                           family = "gaussian",
                           control = myControl)

    test <- predict(hierr_object2,
                    p = 4,
                    pext = 0.05,
                    type = "coefficients",
                    x = xtest,
                    y = ytest,
                    external = ztest,
                    control = myControl)

    expect_equal(drop(predict(hierr_object2,
                              p = 4,
                              pext = 0.05,
                              type = "coefficients",
                              x = xtest,
                              y = ytest,
                              external = ztest,
                              control = myControl)$betas),
                 hierr_object$betas[, 17, 6],
                 tolerance = 1e-6)

    expect_equal(drop(test$beta0), hierr_object$beta0[17, 6], tolerance = 1e-6)
    expect_equal(drop(test$alphas), hierr_object$alphas[ , 17, 6], tolerance = 1e-6)
    expect_equal(drop(test$alpha0), hierr_object$alpha0[17, 6], tolerance = 1e-6)
})

test_that("predict returns right predictions for penalties not already fit by hierr object", {
    myPenalty <- definePenalty(penalty_type = 0,
                               penalty_type_ext = 1,
                               user_penalty = c(2, 1, 0.05),
                               user_penalty_ext = c(0.2, 0.1, 0.05))

    hierr_object <- hierr(x = xtest,
                          y = ytest,
                          external = ztest,
                          family = "gaussian",
                          penalty = myPenalty,
                          control = myControl)

    myPenalty2 <- definePenalty(penalty_type = 0,
                               penalty_type_ext = 1,
                               user_penalty = c(4, 3),
                               user_penalty_ext = c(0.25, 0.15))

    hierr_object2 <- hierr(x = xtest,
                          y = ytest,
                          external = ztest,
                          family = "gaussian",
                          penalty = myPenalty2,
                          control = myControl)

    predy <- cbind(1, xtest) %*% c(hierr_object$beta0[2, 2], hierr_object$betas[, 2, 2])
    expect_equivalent(predict(hierr_object2, p = 1, pext = 0.1, newdata = xtest,
                              x = xtest, y = ytest, external = ztest, penalty = myPenalty, control = myControl),
                      predy,
                      tolerance = 1e-6)

    predy <- cbind(1, xtest) %*% c(hierr_object$beta0[2, 3], hierr_object$betas[, 2, 3])
    expect_equivalent(predict(hierr_object2, p = 1, pext = 0.05, newdata = xtest,
                              x = xtest, y = ytest, external = ztest, penalty = myPenalty, control = myControl),
                      predy,
                      tolerance = 1e-6)

    predy <- cbind(1, xtest) %*% c(hierr_object$beta0[1, 3], hierr_object$betas[, 1, 3])
    expect_equivalent(predict(hierr_object2, p = 2, pext = 0.05, newdata = xtest,
                              x = xtest, y = ytest, external = ztest, penalty = myPenalty, control = myControl),
                      predy,
                      tolerance = 1e-6)
})

test_that("predict returns right predictions for penalties already fit by hierr object, no external data", {
    myPenalty <- definePenalty(penalty_type = 0,
                               penalty_type_ext = 1,
                               user_penalty = c(2, 1, 0.05))

    hierr_object <- hierr(x = xtest,
                          y = ytest,
                          family = "gaussian",
                          intercept = c(T, F),
                          penalty = myPenalty,
                          control = myControl)

    predy1 <- cbind(1, xtest) %*% c(hierr_object$beta0[1, 1], hierr_object$betas[, 1, 1])
    predy2 <- cbind(1, xtest) %*% c(hierr_object$beta0[2, 1], hierr_object$betas[, 2, 1])
    predy3 <- cbind(1, xtest) %*% c(hierr_object$beta0[3, 1], hierr_object$betas[, 3, 1])
    expect_equivalent(predict(hierr_object, p = c(0.05, 1, 2), newdata = xtest), cbind(predy1, predy2, predy3))
})
