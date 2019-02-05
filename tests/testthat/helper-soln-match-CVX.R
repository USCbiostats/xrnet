xtest <- readRDS(file = "testdata/xtest.rds")
ytest <- readRDS(file = "testdata/ytest.rds")
ztest <- readRDS(file = "testdata/ztest.rds")
alphas_cvx_mat <- readRDS(file = "testdata/alphas_cvx_mat.rds")
betas_cvx_mat <- readRDS(file = "testdata/betas_cvx_mat.rds")
sdy <- sqrt(var(ytest) * (length(ytest) - 1) / length(ytest))
ytest_scaled <- ytest / sdy
