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

    if (intr[2]) {
        ptype <- c(penalty$penalty_type, 0, penalty$penalty_type_ext)
        cmult <- c(penalty$custom_multiplier, 0.0, penalty$custom_multiplier_ext)
        lower_cl <- c(control$lower_limits[1:ncol(x)], -Inf, control$lower_limits[(ncol(x) + 1):(ncol(x) + ncol(external))])
        upper_cl <- c(control$upper_limits[1:ncol(x)], Inf, control$upper_limits[(ncol(x) + 1):(ncol(x) + ncol(external))])
    } else {
        ptype <- c(penalty$penalty_type, penalty$penalty_type_ext)
        cmult <- c(penalty$custom_multiplier, penalty$custom_multiplier_ext)
        lower_cl <- control$lower_limits
        upper_cl <- control$upper_limits
    }

    fit <- gaussian_fit(ka = ka,
                        ptype = ptype,
                        nobs = length(y),
                        nvar = ncol(x),
                        nvar_ext = ncol(external),
                        x_ = x,
                        y_ = y,
                        ext_ = external,
                        w = weights,
                        cmult = cmult,
                        lower_cl = lower_cl,
                        upper_cl = upper_cl,
                        ne = control$dfmax[1],
                        ne_ext = control$dfmax[2],
                        nx = control$pmax[1],
                        nx_ext = control$pmax[2],
                        nlam = penalty$num_penalty,
                        nlam_ext = penalty$num_penalty_ext,
                        pratio = penalty$penalty_ratio,
                        pratio_ext = penalty$penalty_ratio_ext,
                        ulam_ = penalty$user_penalty,
                        ulam_ext_ = penalty$user_penalty_ext,
                        thr = control$tolerance,
                        maxit = control$max_iterations,
                        isd = isd[1],
                        isd_ext = isd[2],
                        intr = intr[1],
                        intr_ext = intr[2])
}
