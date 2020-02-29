# linear data
xtest <- readRDS(file = "testdata/xtest.rds")
ytest <- readRDS(file = "testdata/ytest.rds")
ztest <- readRDS(file = "testdata/ztest.rds")

# scaled y
n <- length(ytest)
sd_y <- sqrt(var(ytest) * (n - 1) / n)
ytest_scaled <- ytest / sd_y

# sparse versions of features / external data
xsparse <- Matrix::Matrix(xtest, sparse = TRUE)
zsparse <- Matrix::Matrix(ztest, sparse = TRUE)

# binary data
xtest_binomial <- readRDS(file = "testdata/xtest_binomial.rds")
ytest_binomial <- readRDS(file = "testdata/ytest_binomial.rds")
