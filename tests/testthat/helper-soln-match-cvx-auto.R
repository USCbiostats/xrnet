xtest <- readRDS(file = "testdata/xtest.rds")
ytest <- readRDS(file = "testdata/ytest.rds")
ztest <- readRDS(file = "testdata/ztest.rds")
alphas_cvx_auto <- readRDS(file = "testdata/alphas_cvx_auto.rds")
betas_cvx_auto <- readRDS(file = "testdata/betas_cvx_auto.rds")
n <- length(ytest)
sd_y <- sqrt(var(ytest) * (n - 1) / n)
ytest_scaled <- ytest / sd_y
