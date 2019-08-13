library(parallel)
library(doParallel)
library(bigmemory)
library(glmnet)

context("check computation of CV fold errors")

test_that("gaussian, mse (sequential)", {

    myPenalty <- define_penalty(
        penalty_type = 0,
        num_penalty = 20,
        penalty_type_ext = 1,
        num_penalty_ext = 20
    )

    fit_xrnet <- tune_xrnet(
        x = xtest,
        y = ytest,
        external = ztest,
        family = "gaussian",
        penalty = myPenalty,
        control = list(tolerance = 1e-10),
        loss = "mse",
        foldid = foldid
    )

    expect_equal(
        cv_mean,
        fit_xrnet$cv_mean,
        tolerance = 1e-5,
        check.attribute = FALSE
    )

    fit_xrnet <- tune_xrnet(
        x = xsparse,
        y = ytest,
        external = ztest,
        family = "gaussian",
        penalty = myPenalty,
        control = list(tolerance = 1e-10),
        loss = "mse",
        foldid = foldid
    )

    expect_equal(
        cv_mean,
        fit_xrnet$cv_mean,
        tolerance = 1e-5,
        check.attribute = FALSE
    )

    fit_xrnet <- tune_xrnet(
        x = xtest,
        y = ytest,
        external = zsparse,
        family = "gaussian",
        penalty = myPenalty,
        control = list(tolerance = 1e-10),
        loss = "mse",
        foldid = foldid
    )

    expect_equal(
        cv_mean,
        fit_xrnet$cv_mean,
        tolerance = 1e-5,
        check.attribute = FALSE
    )
})

test_that("gaussian, mse (parallel)", {

    myPenalty <- define_penalty(
        penalty_type = 0,
        num_penalty = 20,
        penalty_type_ext = 1,
        num_penalty_ext = 20
    )

    cl <- makeCluster(2, type = "PSOCK")
    registerDoParallel(cl)

    fit_xrnet <- tune_xrnet(
        x = xtest,
        y = ytest,
        external = ztest,
        family = "gaussian",
        penalty = myPenalty,
        control = list(tolerance = 1e-10),
        loss = "mse",
        foldid = foldid,
        parallel = TRUE
    )

    expect_equal(
        cv_mean,
        fit_xrnet$cv_mean,
        tolerance = 1e-5,
        check.attribute = FALSE
    )

    fit_xrnet <- tune_xrnet(
        x = as.big.matrix(xtest),
        y = ytest,
        external = ztest,
        family = "gaussian",
        penalty = myPenalty,
        control = list(tolerance = 1e-10),
        loss = "mse",
        foldid = foldid,
        parallel = TRUE
    )

    expect_equal(
        cv_mean,
        fit_xrnet$cv_mean,
        tolerance = 1e-5,
        check.attribute = FALSE
    )
})

test_that("gaussian, mae (sequential)", {

    myPenalty <- define_penalty(
        penalty_type = 0,
        num_penalty = 20
    )

    fit_xrnet <- tune_xrnet(
        x = xtest,
        y = ytest,
        family = "gaussian",
        penalty = myPenalty,
        control = list(tolerance = 1e-10),
        loss = "mae",
        foldid = foldid
    )

    fit_glmnet <- cv.glmnet(
        x = xtest,
        y = ytest,
        family = "gaussian",
        alpha = 0,
        foldid = foldid,
        lambda = fit_xrnet$fitted_model$penalty,
        thresh = 1e-10,
        type.measure = "mae"
    )

    expect_equal(unname(fit_glmnet$cvm), drop(fit_xrnet$cv_mean), check.attribute = FALSE)
})

n <- 100
p <- 10
xbin <- matrix(rnorm(n*p), n, p)
b <- rnorm(p)
ybin <- rbinom(n, 1, prob = exp(1 + xbin %*% b) / (1 + exp(1 + xbin %*% b)))
foldid_bin <- sample(rep(seq(5), length = n))

test_that("binomial, auc (sequential)", {

    myPenalty <- define_penalty(
        penalty_type = 0,
        num_penalty = 20,
        penalty_ratio = 0.001
    )

    fit_xrnet <- tune_xrnet(
        x = xbin,
        y = ybin,
        family = "binomial",
        penalty = myPenalty,
        control = list(tolerance = 1e-10),
        loss = "auc",
        foldid = foldid_bin
    )

    fit_glmnet <- cv.glmnet(
        x = xbin,
        y = ybin,
        family = "binomial",
        alpha = 0,
        foldid = foldid_bin,
        lambda = fit_xrnet$fitted_model$penalty,
        thresh = 1e-10,
        type.measure = "auc"
    )

    expect_equal(unname(fit_glmnet$cvm), drop(fit_xrnet$cv_mean), check.attribute = FALSE)
})

test_that("binomial, deviance (sequential)", {

    myPenalty <- define_penalty(
        penalty_type = 0,
        num_penalty = 20,
        penalty_ratio = 0.001
    )

    fit_xrnet <- tune_xrnet(
        x = xbin,
        y = ybin,
        family = "binomial",
        penalty = myPenalty,
        control = list(tolerance = 1e-10),
        loss = "deviance",
        foldid = foldid_bin
    )

    fit_glmnet <- cv.glmnet(
        x = xbin,
        y = ybin,
        family = "binomial",
        alpha = 0,
        foldid = foldid_bin,
        lambda = fit_xrnet$fitted_model$penalty,
        thresh = 1e-10,
        type.measure = "deviance"
    )

    expect_equal(unname(fit_glmnet$cvm), drop(fit_xrnet$cv_mean), check.attribute = FALSE)
})
