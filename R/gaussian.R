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

    fit <- gaussian_fit(x_ = x,
                        y_ = y,
                        ext_ = external,
                        nobs = nrow(x),
                        nvar = ncol(x),
                        nvar_ext = ncol(external),
                        w = weights,
                        ptype = penalty$ptype,
                        cmult = penalty$cmult,
                        lower_cl = control$lower_limits,
                        upper_cl = control$upper_limits,
                        ne = control$dfmax,
                        nx = control$pmax,
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
