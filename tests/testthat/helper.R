xtest <- readRDS(file = "testdata/xtest.rds")
ytest <- readRDS(file = "testdata/ytest.rds")
ztest <- readRDS(file = "testdata/ztest.rds")
n <- length(ytest)
sd_y <- sqrt(var(ytest) * (n - 1) / n)
ytest_scaled <- ytest / sd_y
xsparse <- Matrix::Matrix(xtest, sparse = TRUE)
zsparse <- Matrix::Matrix(ztest, sparse = TRUE)
