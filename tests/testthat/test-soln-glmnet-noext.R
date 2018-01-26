context("compare coefficients to glmnet when no external data")

##### Code used to generate data files for all tests #####

# set.seed(123)
# n <- 100
# p <- 50
# q <- 5

# meanx_true <- rep(0, p)
# covx_true <- matrix(NA, nrow = p, ncol = p)
# for (i in 1:p) {
#     for (j in 1:p) {
#         covx_true[i, j] <- 0.5^abs(i-j)
#     }
# }

# meanz_true <- rep(0, q)
# covz_true <- matrix(NA, nrow = q, ncol = q)
# for (i in 1:q) {
#     for (j in 1:q) {
#         covz_true[i, j] <- 0.5^abs(i-j)
#     }
# }

# a0 <- 0.01
# a <- c(0.1, -0.1, rep(0, q - 2))
# z <- mvrnorm(n = p, mu = meanz_true, Sigma = covz_true)
# b <- drop(a0 + z %*% a + 0.2*rnorm(p))
# x <- mvrnorm(n = n, mu = meanx_true, Sigma = covx_true)
# y <- drop(x %*% b + rnorm(n))

# mean_x <- colMeans(x)
# var_x <- apply(x, 2, var) * (n-1) / n
# sd_x <- sqrt(var_x)

# xscaled <- matrix(NA, nrow = n, ncol = p)
# for (i in 1:p) {
#     xscaled[, i] <- (x[, i] - mean_x[i]) / sd_x[i]
# }

# mean_z <- colMeans(z)
# sd_z <- sqrt(apply(z, 2, var) * (p-1) / p)

# zscaled <- matrix(NA, nrow = p, ncol = q)
# for (i in 1:q) {
#     zscaled[, i] <- (z[, i] - mean_z[i]) / sd_z[i]
# }

load(file = "Test-Data/x.Rdata")
load(file = "Test-Data/y.Rdata")
load(file = "Test-Data/z.Rdata")

myPenalty <- definePenalty(penalty_type = 0, num_penalty = 100)
myControl <- list(tolerance = 1e-20)

# Ridge Regression #

test_that("x standardized, intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", alpha = 0, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), penalty = myPenalty, control = myControl)$betas[, , 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), penalty = myPenalty, control = myControl)$beta0[,1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, alpha = 0, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$betas[, , 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$beta0[,1],
                 tolerance = 1e-5)
})

test_that("x standardized, NO intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", intercept = FALSE, alpha = 0, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), penalty = myPenalty, control = myControl)$betas[, , 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), penalty = myPenalty, control = myControl)$beta0[,1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, NO intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, intercept = FALSE, alpha = 0, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$betas[, , 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$beta0[,1],
                 tolerance = 1e-5)
})

# Lasso Regression #

myPenalty <- definePenalty(penalty_type = 1, num_penalty = 100)

test_that("x standardized, intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", alpha = 1, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), penalty = myPenalty, control = myControl)$betas[, 1:80, 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), penalty = myPenalty, control = myControl)$beta0[1:80, 1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, alpha = 1, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$betas[, 1:81, 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$beta0[1:81, 1],
                 tolerance = 1e-5)
})

test_that("x standardized, NO intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", intercept = FALSE, alpha = 1, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), penalty = myPenalty, control = myControl)$betas[, 1:80, 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), penalty = myPenalty, control = myControl)$beta0[1:80, 1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, NO intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, intercept = FALSE, alpha = 1, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$betas[, 1:81, 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$beta0[1:81, 1],
                 tolerance = 1e-5)
})

# Elastic Net Regression #

myPenalty <- definePenalty(penalty_type = 0.5, num_penalty = 100)

test_that("x standardized, intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", alpha = 0.5, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), penalty = myPenalty, control = myControl)$betas[, 1:81, 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), penalty = myPenalty, control = myControl)$beta0[1:81, 1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, alpha = 0.5, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$betas[, 1:82, 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(T, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$beta0[1:82, 1],
                 tolerance = 1e-5)
})

test_that("x standardized, NO intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", intercept = FALSE, alpha = 0.5, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), penalty = myPenalty, control = myControl)$betas[, 1:82, 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), penalty = myPenalty, control = myControl)$beta0[1:82, 1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, NO intercept",{

    coef_glmnet <- glmnet::glmnet(x = x, y = y, family = "gaussian", standardize = FALSE, intercept = FALSE, alpha = 0.5, thresh = 1e-20)
    betas_glmnet <- unname(as.matrix(coef_glmnet$beta))
    b0_glmnet <- unname(coef_glmnet$a0)

    expect_equal(betas_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$betas[, 1:83, 1],
                 tolerance = 1e-5)

    expect_equal(b0_glmnet,
                 hierr(x = x, y = y, family = "gaussian", intercept = c(F, F), standardize = c(F, T), penalty = myPenalty, control = myControl)$beta0[1:83, 1],
                 tolerance = 1e-5)
})
