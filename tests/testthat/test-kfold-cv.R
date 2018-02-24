context("test k-fold cross-validation function")
load(file = "Test-Data/x.Rdata")
load(file = "Test-Data/y.Rdata")
load(file = "Test-Data/z.Rdata")
sd_y <- sqrt(var(y) * (length(y)-1) / length(y))

test_that("obtain correct min(penalty) compared to glmnet (no external data) -- Ridge",{
    # glmnet code used to find min(lambda)
    #set.seed(123)
    #cv_glmnet <- cv.glmnet(x, y, family = "gaussian", nfolds = 5, alpha = 0, thresh = 1e-15, keep = T)
    set.seed(123)
    myPenalty <- definePenalty(penalty_type = 0, num_penalty = 100)
    myControl <- list(tolerance = 1e-15)
    expect_equal(cvhierr(x = x, y = y, family = "gaussian", intercept = c(T, F), penalty = myPenalty, control = myControl)$minl1, 0.2445511, tolerance = 1e-6)

})

test_that("obtain correct min(penalty) compared to glmnet (no external data) -- Lasso",{
    # glmnet code used to find min(lambda)
    # set.seed(123)
    # cv_glmnet <- cv.glmnet(x, y, family = "gaussian", nfolds = 5, alpha = 1, thresh = 1e-15)
    set.seed(123)
    myPenalty <- definePenalty(penalty_type = 1, num_penalty = 100)
    myControl <- list(tolerance = 1e-15)
    expect_equal(cvhierr(x = x, y = y, family = "gaussian", intercept = c(T, F), penalty = myPenalty, control = myControl)$minl1, 0.0371696, tolerance = 1e-6)

})

test_that("obtain correct min(penalty) compared to glmnet (no external data) -- Elastic Net",{
    # glmnet code used to find min(lambda)
    #set.seed(123)
    #cv_glmnet <- cv.glmnet(x, y, family = "gaussian", nfolds = 5, alpha = 0.5, thresh = 1e-15)
    set.seed(123)
    myPenalty <- definePenalty(penalty_type = 0.5, num_penalty = 100)
    myControl <- list(tolerance = 1e-15)
    expect_equal(cvhierr(x = x, y = y, family = "gaussian", intercept = c(T, F), penalty = myPenalty, control = myControl)$minl1, 0.06773511, tolerance = 1e-6)

})
