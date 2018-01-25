# Function to prepare Gaussian outcome for model fitting

gaussian <- function(x,
                     nr_x,
                     nc_x,
                     y,
                     external,
                     nc_ext,
                     penalty,
                     ptype,
                     cmult,
                     weights,
                     isd,
                     intr,
                     control) {

    y <- as.double(drop(y))

    fit <- gaussian_fit(ptype = penalty$ptype,
                        nobs = nr_x,
                        nvar = nc_x,
                        nvar_ext = nc_ext,
                        x_ = x,
                        y_ = y,
                        ext_ = external,
                        w = weights,
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
