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
        thresh = 1e-15,
        alpha = 0,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 0, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0), tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas), tolerance = 1e-5)
})

test_that("x NOT standardized, intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-15,
        alpha = 0,
        standardize = FALSE,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 0, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        standardize = c(F, F),
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0), tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas), tolerance = 1e-5)
})

test_that("x standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-15,
        alpha = 0,
        intercept = FALSE,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 0, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0), tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas), tolerance = 1e-5)
})

test_that("x NOT standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-15,
        alpha = 0,
        intercept = FALSE,
        standardize = FALSE,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 0, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        standardize = c(F,F),
        control = xrnet.control(tolerance = 1e-15)
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
        thresh = 1e-15,
        alpha = 1,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 1, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x NOT standardized, intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-15,
        alpha = 1,
        standardize = FALSE,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 1, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        standardize = c(F, F),
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-15,
        alpha = 1,
        intercept = FALSE,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 1, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x NOT standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-15,
        alpha = 1,
        intercept = FALSE,
        standardize = FALSE,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 1, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        standardize = c(F,F),
        control = xrnet.control(tolerance = 1e-15)
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
        thresh = 1e-15,
        alpha = 0.5,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 0.5, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x NOT standardized, intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-15,
        alpha = 0.5,
        standardize = FALSE,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 0.5, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        standardize = c(F, F),
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-15,
        alpha = 0.5,
        intercept = FALSE,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 0.5, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

test_that("x NOT standardized, NO intercept",{

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-15,
        alpha = 0.5,
        intercept = FALSE,
        standardize = FALSE,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(penalty_type = 0.5, num_penalty = 100, penalty_ratio = 0.01)

    fit_xrnet <- xrnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        penalty = myPenalty,
        intercept = c(F, F),
        standardize = c(F,F),
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(unname(fit_glmnet$a0), drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)], tolerance = 1e-5)
    expect_equal(unname(as.matrix(fit_glmnet$beta)), drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)], tolerance = 1e-5)
})

# Elastic Net - No Penalty on 1st two variables #

test_that("x NOT standardized, intercept",{

    pf <- c(rep(0, 2), rep(1, NCOL(xtest) - 2))

    fit_glmnet <- glmnet(
        x = xtest,
        y = ytest,
        family = "binomial",
        thresh = 1e-15,
        alpha = 0.5,
        penalty.factor = pf,
        lambda.min.ratio = 0.01
    )

    myPenalty <- define_penalty(
        penalty_type = 0.5,
        user_penalty = fit_glmnet$lambda,
        custom_multiplier = (NCOL(xtest) / sum(pf)) * rep(1, NCOL(xtest) - 2)
    )

    fit_xrnet <- xrnet(
        x = xtest[, -c(1, 2)],
        y = ytest,
        unpen = xtest[, c(1, 2)],
        family = "binomial",
        penalty = myPenalty,
        control = xrnet.control(tolerance = 1e-15)
    )

    expect_equal(
        unname(fit_glmnet$a0),
        drop(fit_xrnet$beta0)[1:length(fit_glmnet$a0)],
        tolerance = 1e-5
    )

    expect_equal(
        unname(as.matrix(fit_glmnet$beta)[1:2, ]),
        drop(fit_xrnet$gammas)[, 1:length(fit_glmnet$a0)],
        tolerance = 1e-5
    )

    expect_equal(
        unname(as.matrix(fit_glmnet$beta)[3:10,]),
        drop(fit_xrnet$betas)[, 1:length(fit_glmnet$a0)],
        tolerance = 1e-5
    )
})

