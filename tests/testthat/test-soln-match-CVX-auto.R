context("compare coefficent estimates to CVX (automatic penalty)")

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

myPenalty <- hierr::definePenalty(0, 1)

test_that("x and ext standardized, both intercepts",{

    load(file = "Test-Data/alphas_cvx13.Rdata")
    load(file = "Test-Data/betas_cvx13.Rdata")

    expect_equal(alphas_cvx13,
                 hierr(x = x, y = y, ext = z, family = "gaussian", penalty = myPenalty, control = list(tolerance = 1e-20))$alphas[1:5, 4, 2],
                 tolerance = 1e-5)
    expect_equal(betas_cvx13,
                 hierr(x = x, y = y, ext = z, family = "gaussian", penalty = myPenalty, control = list(tolerance = 1e-20))$betas[1:50, 4, 2],
                 tolerance = 1e-5)
})
