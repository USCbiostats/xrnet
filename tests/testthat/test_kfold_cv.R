context("test k-fold cross-validation function")
library(hierr)

xtest <- readRDS(file = "testdata/xtest.rds")
ytest <- readRDS(file = "testdata/ytest.rds")
ztest <- readRDS(file = "testdata/ztest.rds")

myControl <- list(tolerance = 1e-15, earlyStop = FALSE)

test_that("obtain correct min(penalty) compared to glmnet (no external data) -- Ridge",{
    # glmnet code used to find min(lambda)
    #set.seed(123)
    #cv_glmnet <- cv.glmnet(x, y, family = "gaussian", nfolds = 5, alpha = 0, thresh = 1e-15, keep = T)
    set.seed(123)
    myPenalty <- definePenalty(penalty_type = 0, num_penalty = 100)
    expect_equal(cvhierr(x = xtest,
                         y = ytest,
                         family = "gaussian",
                         intercept = c(T, F),
                         penalty = myPenalty,
                         control = myControl)$minl1,
                 0.5597391,
                 tolerance = 1e-5)

})

test_that("obtain correct min(penalty) compared to glmnet (no external data) -- Lasso",{
    # glmnet code used to find min(lambda)
    #set.seed(123)
    #cv_glmnet <- cv.glmnet(x, y, family = "gaussian", nfolds = 5, alpha = 1, thresh = 1e-15)
    set.seed(123)
    myPenalty <- definePenalty(penalty_type = 1, num_penalty = 100)
    expect_equal(cvhierr(x = xtest,
                         y = ytest,
                         family = "gaussian",
                         intercept = c(T, F),
                         penalty = myPenalty,
                         control = myControl)$minl1,
                 0.08507537,
                 tolerance = 1e-5)

})

test_that("obtain correct min(penalty) compared to glmnet (no external data) -- Elastic Net",{
    # glmnet code used to find min(lambda)
    #set.seed(123)
    #cv_glmnet <- cv.glmnet(x, y, family = "gaussian", nfolds = 5, alpha = 0.5, thresh = 1e-15)
    set.seed(123)
    myPenalty <- definePenalty(penalty_type = 0.5, num_penalty = 100)
    expect_equal(cvhierr(x = xtest,
                         y = ytest,
                         family = "gaussian",
                         intercept = c(T, F),
                         penalty = myPenalty,
                         control = myControl)$minl1,
                 0.155035,
                 tolerance = 1e-5)

})
