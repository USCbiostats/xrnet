# Function to prepare Gaussian outcome for model fitting

gaussian <- function(x,
                     y,
                     external,
                     penalty,
                     weights,
                     isd,
                     intr,
                     control) {

    y <- as.double(drop(y))
    ka <- as.integer(1)

    #fit <- gaussian_fit(ka = ka,
    #                    ptype = c(penalty$penalty_type, penalty$penalty_type_ext),
    #                    nobs = length(y),
    #                    nvar = ncol(x),
    #                    nvar_ext = ncol(external),
    #                    w = weights,
    #                    cmult = c(penalty$custom_multiplier, 1, penalty$custom_multiplier_ext),
    #                    lower_cl = control$lower_limits,
    #                    upper_cl = control$upper_limits,
    #                    ne = control$dfmax,
    #                    nx = control$pmax,
    #                    nlam = penalty$num_penalty,
    #                    nlam_ext = penalty$num_penalty_ext,
    #                    flmin = penalty$penalty_ratio,
    #                    flmin_ext = penalty$penalty_ratio_ext,
    #                    ulam_ = penalty$user_penalty,
    #                    ulam_ext_ = penalty$user_penalty_ext,
    #                    thr = control$tolerance,
    #                    maxit = control$max_iterations,
    #                    isd = isd[1],
    #                    isd_ext = isd[2],
    #                    intr = intr[1],
    #                    intr_ext = intr[2])
}
