library(glmnet)

context("compare coefficients to glmnet when no external data (binomial)")

n <- 100
p <- 10
xtest <- matrix(rnorm(n*p), n, p)
b <- rnorm(p)
ytest <- rbinom(n, 1, prob = exp(1 + xtest %*% b) / (1 + exp(1 + xtest %*% b)))

# Ridge Regression #

test_that("x standardized, intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 0
    )

    myPenalty <- define_penalty(penalty_type = 0, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0), tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas), tolerance = 1e-5)
})

test_that("x NOT standardized, intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 0,
        standardize = FALSE
    )

    myPenalty <- define_penalty(penalty_type = 0, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        standardize = c(F, F),
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0), tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas), tolerance = 1e-5)
})

test_that("x standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 0,
        intercept = FALSE
    )

    myPenalty <- define_penalty(penalty_type = 0, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0), tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas), tolerance = 1e-5)
})

test_that("x NOT standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 0,
        intercept = FALSE,
        standardize = FALSE
    )

    myPenalty <- define_penalty(penalty_type = 0, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        standardize = c(F,F),
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0), tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas), tolerance = 1e-5)
})

# Lasso Regression #

test_that("x standardized, intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 1
    )

    myPenalty <- define_penalty(penalty_type = 1, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x NOT standardized, intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 1,
        standardize = FALSE
    )

    myPenalty <- define_penalty(penalty_type = 1, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        standardize = c(F, F),
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 1,
        intercept = FALSE
    )

    myPenalty <- define_penalty(penalty_type = 1, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x NOT standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 1,
        intercept = FALSE,
        standardize = FALSE
    )

    myPenalty <- define_penalty(penalty_type = 1, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        standardize = c(F,F),
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

# Elastic Net Regression #

test_that("x standardized, intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 0.5
    )

    myPenalty <- define_penalty(penalty_type = 0.5, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x NOT standardized, intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 0.5,
        standardize = FALSE
    )

    myPenalty <- define_penalty(penalty_type = 0.5, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        standardize = c(F, F),
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 0.5,
        intercept = FALSE
    )

    myPenalty <- define_penalty(penalty_type = 0.5, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x NOT standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-12,
        alpha = 0.5,
        intercept = FALSE,
        standardize = FALSE
    )

    myPenalty <- define_penalty(penalty_type = 0.5, num_penalty = 100)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        standardize = c(F,F),
        control = xrnet.control(tolerance = 1e-12)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})
