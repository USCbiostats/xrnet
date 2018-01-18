context("test error function generator for CV function")

load("Test-Data/x.Rdata")
load("Test-Data/y.Rdata")
load("Test-Data/betas_cvx.Rdata")
load("Test-Data/betas_cvx2.Rdata")
betas <- cbind(betas_cvx, betas_cvx2)
w <- rep(1.0, length(y))

test_that("mse function correctly generated",{
    expect_equal(hierr:::error_match("gaussian", "mse")(betas, y, x, w)[ ,1], drop(w * (y - x %*% betas_cvx)^2))
    expect_equal(hierr:::error_match("gaussian", "mse")(betas, y, x, w)[ ,2], drop(w * (y - x %*% betas_cvx2)^2))
})

test_that("deviance function correctly generated",{
    expect_equal(hierr:::error_match("gaussian", "deviance")(betas, y, x, w)[ ,1], drop(w * (y - x %*% betas_cvx)^2))
    expect_equal(hierr:::error_match("gaussian", "deviance")(betas, y, x, w)[ ,2], drop(w * (y - x %*% betas_cvx2)^2))
})

test_that("mae function correctly generated",{
    expect_equal(hierr:::error_match("gaussian", "mae")(betas, y, x, w)[ ,1], drop(w * abs(y - x %*% betas_cvx)))
    expect_equal(hierr:::error_match("gaussian", "mae")(betas, y, x, w)[ ,2], drop(w * abs(y - x %*% betas_cvx2)))
})
