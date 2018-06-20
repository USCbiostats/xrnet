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

xtest <- readRDS(file = "testdata/xtest.rds")
ytest <- readRDS(file = "testdata/ytest.rds")
ztest <- readRDS(file = "testdata/ztest.rds")
alphas_cvx_auto <- readRDS(file = "testdata/alphas_cvx_auto.rds")
betas_cvx_auto <- readRDS(file = "testdata/betas_cvx_auto.rds")

myPenalty <- definePenalty(penalty_type = 0, penalty_type_ext = 1)
myControl <- list(tolerance = 1e-15, earlyStop = TRUE)

test_that("x and ext standardized, both intercepts, auto",{

    expect_equal(alphas_cvx_auto[, 1],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 1],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext is standardized, both intercepts, auto",{

    expect_equal(alphas_cvx_auto[, 2],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, T),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 2],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, T),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})

test_that("x is standardized, ext NOT standardized, both intercepts, auto",{

    expect_equal(alphas_cvx_auto[, 3],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(T, F),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 3],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(T, F),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext NOT standardized, both intercepts, auto",{

    expect_equal(alphas_cvx_auto[, 4],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, F),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 4],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, F),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})


############################# No 2nd level intercept ########################################

test_that("x standardized, ext standardized, no 2nd level, auto",{

    expect_equal(alphas_cvx_auto[, 5],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(T, T),
                       intercept = c(T, F),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 5],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(T, T),
                       intercept = c(T, F),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext is standardized, no 2nd level, auto",{

    expect_equal(alphas_cvx_auto[, 6],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, T),
                       intercept = c(T, F),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 6],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, T),
                       intercept = c(T, F),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})

test_that("x is standardized, ext NOT standardized, no 2nd level, auto",{

    expect_equal(alphas_cvx_auto[, 7],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(T, F),
                       intercept = c(T, F),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 7],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(T, F),
                       intercept = c(T, F),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext NOT standardized, no 2nd level, auto",{

    expect_equal(alphas_cvx_auto[, 8],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, F),
                       intercept = c(T, F),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 8],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, F),
                       intercept = c(T, F),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})

############################# No 2nd level intercept ########################################

test_that("x standardized, ext standardized, no 1st level, auto",{

    expect_equal(alphas_cvx_auto[, 9],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(T, T),
                       intercept = c(F, T),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 9],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(T, T),
                       intercept = c(F, T),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext is standardized, no 1st level, auto",{

    expect_equal(alphas_cvx_auto[, 10],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, T),
                       intercept = c(F, T),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 10],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, T),
                       intercept = c(F, T),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})

test_that("x is standardized, ext NOT standardized, no 1st level, auto",{

    expect_equal(alphas_cvx_auto[, 11],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(T, F),
                       intercept = c(F, T),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 11],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(T, F),
                       intercept = c(F, T),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})

test_that("x NOT standardized, ext NOT standardized, no 1st level, auto",{

    expect_equal(alphas_cvx_auto[, 12],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, F),
                       intercept = c(F, T),
                       penalty = myPenalty,
                       control = myControl)$alphas[, 4, 2],
                 tolerance = 1e-5)

    expect_equal(betas_cvx_auto[, 12],
                 hierr(x = xtest,
                       y = ytest,
                       external = ztest,
                       family = "gaussian",
                       standardize = c(F, F),
                       intercept = c(F, T),
                       penalty = myPenalty,
                       control = myControl)$betas[, 4, 2],
                 tolerance = 1e-5)
})
