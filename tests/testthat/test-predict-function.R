context("test predict function works correctly")

load(file = "Test-Data/x.Rdata")
load(file = "Test-Data/y.Rdata")
load(file = "Test-Data/z.Rdata")

test_that("predict returns estimates for penalties already fit by hierr object", {
    hierr_object <- hierr(x, y, z, family = "gaussian", penalty = definePenalty(0, 1, user_penalty = c(2, 1, 0.05), user_penalty_ext = c(0.2, 0.1, 0.05)))
    test <- predict(hierr_object, p = 1, pext = 0.05, type = "coefficients")
    expect_identical(drop(predict(hierr_object, p = 1, pext = 0.05, type = "coefficients")$betas), hierr_object$betas[, 2, 3])
    expect_identical(drop(test$beta0), hierr_object$beta0[2, 3])
    expect_identical(drop(test$alphas), hierr_object$alphas[ , 2, 3])
    expect_identical(drop(test$alpha0), hierr_object$alpha0[2, 3])
})


test_that("predict returns estimates for penalties notfit by hierr object", {
    hierr_object <- hierr(x, y, z, family = "gaussian", penalty = definePenalty(0, 1, user_penalty = 1:20, user_penalty_ext = seq(0.01, 0.1, 0.01)), control = list(tolerance = 1e-20))
    hierr_object2 <- hierr(x, y, z, family = "gaussian", control = list(tolerance = 1e-20))
    test <- predict(hierr_object2, p = 4, pext = 0.05, type = "coefficients", x = x, y = y, external = z)
    expect_equal(drop(predict(hierr_object2, p = 4, pext = 0.05, type = "coefficients", x = x, y = y, external = z)$betas), hierr_object$betas[, 17, 6], tolerance = 1e-6)
    expect_equal(drop(test$beta0), hierr_object$beta0[17, 6], tolerance = 1e-6)
    expect_equal(drop(test$alphas), hierr_object$alphas[ , 17, 6], tolerance = 1e-6)
    expect_equal(drop(test$alpha0), hierr_object$alpha0[17, 6], tolerance = 1e-6)
})
