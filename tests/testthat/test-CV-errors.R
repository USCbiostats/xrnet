context("check computation of CV fold errors")

test_that("gaussian", {

    myPenalty <- define_penalty(penalty_type = 0,
                                num_penalty = 20,
                                penalty_type_ext = 1,
                                num_penalty_ext = 20)

    expect_equal(
        cv_mean,
        cv_hierr(x = xtest,
                 y = ytest,
                 external = ztest,
                 family = "gaussian",
                 penalty = myPenalty,
                 control = list(tolerance = 1e-10),
                 loss = "mse",
                 foldid = foldid)$cv_mean,
        tolerance = 1e-5,
        check.attribute = FALSE
    )
})
