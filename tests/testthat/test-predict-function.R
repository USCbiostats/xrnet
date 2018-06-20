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
    expect_identical(drop(predict(hierr_object,
                                  p = 1,
                                  pext = 0.05,
                                  type = "coefficients")$betas),
                     hierr_object$betas[, 2, 3])
    expect_identical(drop(test$beta0), hierr_object$beta0[2, 3])
    expect_identical(drop(test$alphas), hierr_object$alphas[ , 2, 3])
    expect_identical(drop(test$alpha0), hierr_object$alpha0[2, 3])
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
