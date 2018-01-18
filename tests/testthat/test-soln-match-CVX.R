context("compare coefficent estimates to CVX")

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

myPenalty <- hierr::definePenalty(0, 1, user_penalty = 1, user_penalty_ext = 0.1)

test_that("x and ext standardized, both intercepts",{

    load(file = "Test-Data/alphas_cvx.Rdata")
    load(file = "Test-Data/betas_cvx.Rdata")

    expect_equal(alphas_cvx,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)
    expect_equal(betas_cvx,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext standardized, both intercepts",{

    load(file = "Test-Data/alphas_cvx2.Rdata")
    load(file = "Test-Data/betas_cvx2.Rdata")

    expect_equal(alphas_cvx2,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, standardize = c(FALSE, TRUE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx2,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, standardize = c(FALSE, TRUE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x standardized, ext NOT standardized, both intercepts",{

    load(file = "Test-Data/alphas_cvx3.Rdata")
    load(file = "Test-Data/betas_cvx3.Rdata")

    expect_equal(alphas_cvx3,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, standardize = c(TRUE, FALSE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx3,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, standardize = c(TRUE, FALSE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext NOT standardized, both intercepts",{

    load(file = "Test-Data/alphas_cvx4.Rdata")
    load(file = "Test-Data/betas_cvx4.Rdata")

    expect_equal(alphas_cvx4,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, standardize = c(FALSE, FALSE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx4,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, standardize = c(FALSE, FALSE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x standardized, ext standardized, no 2nd level intercept",{

    load(file = "Test-Data/alphas_cvx5.Rdata")
    load(file = "Test-Data/betas_cvx5.Rdata")

    expect_equal(alphas_cvx5,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(TRUE, FALSE), standardize = c(TRUE, TRUE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx5,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(TRUE, FALSE), standardize = c(TRUE, TRUE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext standardized, no 2nd level intercept",{

    load(file = "Test-Data/alphas_cvx6.Rdata")
    load(file = "Test-Data/betas_cvx6.Rdata")

    expect_equal(alphas_cvx6,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(TRUE, FALSE), standardize = c(FALSE, TRUE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx6,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(TRUE, FALSE), standardize = c(FALSE, TRUE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x standardized, ext NOT standardized, no 2nd level intercept",{

    load(file = "Test-Data/alphas_cvx7.Rdata")
    load(file = "Test-Data/betas_cvx7.Rdata")

    expect_equal(alphas_cvx7,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(TRUE, FALSE), standardize = c(TRUE, FALSE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx7,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(TRUE, FALSE), standardize = c(TRUE, FALSE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext NOT standardized, no 2nd level intercept",{

    load(file = "Test-Data/alphas_cvx8.Rdata")
    load(file = "Test-Data/betas_cvx8.Rdata")

    expect_equal(alphas_cvx8,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(TRUE, FALSE), standardize = c(FALSE, FALSE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx8,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(TRUE, FALSE), standardize = c(FALSE, FALSE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x standardized, ext standardized, no 1st level intercept",{

    load(file = "Test-Data/alphas_cvx9.Rdata")
    load(file = "Test-Data/betas_cvx9.Rdata")

    expect_equal(alphas_cvx9,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(FALSE, TRUE), standardize = c(TRUE, TRUE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx9,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(FALSE, TRUE), standardize = c(TRUE, TRUE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext standardized, no 1st level intercept",{

    load(file = "Test-Data/alphas_cvx10.Rdata")
    load(file = "Test-Data/betas_cvx10.Rdata")

    expect_equal(alphas_cvx10,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(FALSE, TRUE), standardize = c(FALSE, TRUE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx10,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(FALSE, TRUE), standardize = c(FALSE, TRUE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x standardized, ext NOT standardized, no 1st level intercept",{

    load(file = "Test-Data/alphas_cvx11.Rdata")
    load(file = "Test-Data/betas_cvx11.Rdata")

    expect_equal(alphas_cvx11,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(FALSE, TRUE), standardize = c(TRUE, FALSE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx11,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(FALSE, TRUE), standardize = c(TRUE, FALSE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext NOT standardized, no 1st level intercept",{

    load(file = "Test-Data/alphas_cvx12.Rdata")
    load(file = "Test-Data/betas_cvx12.Rdata")

    expect_equal(alphas_cvx12,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(FALSE, TRUE), standardize = c(FALSE, FALSE), control = list(tolerance = 1e-20))$alphas[1:5, 1, 1],
                 tolerance = 1e-5)

    expect_equal(betas_cvx12,
                 hierr(x = x, y = y, ext = z, penalty = myPenalty, intercept = c(FALSE, TRUE), standardize = c(FALSE, FALSE), control = list(tolerance = 1e-20))$betas[1:50, 1, 1],
                 tolerance = 1e-5)
})
