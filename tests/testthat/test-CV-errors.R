context("check computation of CV fold errors")

test_that("gaussian", {

    myPenalty <- define_penalty(penalty_type = 0,
                                num_penalty = 20,
                                penalty_type_ext = 1,
                                num_penalty_ext = 20)

    expect_equal(
        cv_mean,
        tune_xrnet(x = xtest,
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

#myPenalty <- define_penalty(penalty_type = 0,
#                            penalty_type_ext = 1,
#                            user_penalty = check$fitted_model$penalty,
#                            user_penalty_ext = check$fitted_model$penalty_ext)

#errormat <- matrix(NA, nrow = 400, ncol = 5)
#for (k in 1:5) {
#
#    fit_fold <- xrnet(x = xtest[foldid != k,],
#                      y = ytest[foldid != k],
#                      external = ztest,
#                      family = "gaussian",
#                      penalty = myPenalty,
#                      control = list(tolerance = 1e-10))
#
#    betas <- rbind(as.vector(t(fit_fold$beta0)),
#                   `dim<-`(aperm(fit_fold$betas, c(1, 3, 2)),
#                           c(dim(fit_fold$betas)[1],
#                             dim(fit_fold$betas)[2] * dim(fit_fold$betas)[3])))
#
#    predy <- cbind(1, xtest[foldid == k, ]) %*% betas
#
#    errormat[, k] <- apply(predy, 2, function(p) mean((ytest[foldid == k] - p)^2))
#}
#cv_mean <- rowMeans(errormat)
#cv_mean <- matrix(cv_mean, nrow = 20, byrow = TRUE)
#all.equal(check$cv_mean, cv_mean, check.attributes = F)



