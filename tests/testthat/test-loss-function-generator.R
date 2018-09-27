context("test error function generator for CV function")

xtest <- readRDS(file = "testdata/xtest.rds")
ytest <- readRDS(file = "testdata/ytest.rds")
betas_cvx_mat <- readRDS(file = "testdata/betas_cvx_mat.rds")

betas <- betas_cvx_mat[,c(1:2)]
w <- rep(1.0, length(ytest))

test_that("mse function correctly generated",{
    expect_equal(hierr:::error_match("gaussian", "mse")(betas, ytest, xtest, w),
                 colSums((w / sum(w) * (ytest - xtest %*% betas)^2)))
})

test_that("deviance function correctly generated",{
    expect_equal(hierr:::error_match("gaussian", "deviance")(betas, ytest, xtest, w),
                 colSums((w / sum(w) * (ytest - xtest %*% betas)^2)))
})

test_that("mae function correctly generated",{
    expect_equal(hierr:::error_match("gaussian", "mae")(betas, ytest, xtest, w),
                 colSums((w / sum(w) * abs(ytest - xtest %*% betas))))
})
