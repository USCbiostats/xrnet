#' @useDynLib hierr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats predict
NULL

#' Fit hierarchical regularized regression model
#'
#' @description Fits hierarchical regularized regression model that enables the incorporation of external data
#' for the predictor variables. Both the predictor variables and external data can be regularized
#' by the most common penalties (lasso, ridge, elastic net) and we have included an addtional "quantile" penalty.
#' Solutions are computed across a two-dimensional grid of penalties (a separate penalty path is computed
#' for the predictors and external variables). Currently support linear regression, future extensions to
#' a GLM framework will be implemented in the next major update.
#'
#' @param x predictor design matrix of dimension \eqn{n x p}, matrix options include:
#' \itemize{
#'    \item big.matrix
#'    \item filebacked.big.matrix
#'    \item sparse matrix (dgCMatrix)
#' }
#' @param y outcome vector of length \eqn{n}
#' @param external (optional) external data design matrix of dimension \eqn{p x q}, matrix options include:
#' \itemize{
#'     \item matrix
#'     \item sparse matrix (dgCMatrix)
#' }
#' @param unpen (optional) unpenalized predictor design matrix
#' @param family error distribution for outcome variable, options include:
#' \itemize{
#'     \item gaussian
#'     \item binomial
#' }
#' @param penalty specifies regularization object for x and external. See \code{\link{definePenalty}} for more details.
#' @param weights optional vector of observation-specific weights. Default is 1 for all observations.
#' @param standardize indicates whether x and/or external should be standardized. Default is c(TRUE, TRUE).
#' @param intercept indicates whether an intercept term is included for x and/or external. Default is c(TRUE, FALSE).
#' @param control specifies hierr control object. See \code{\link{hierr.control}} for more details.
#' @return A list of class \code{hierr} with components
#' \item{beta0}{matrix of first-level intercepts indexed by penalty values}
#' \item{betas}{3-dimensional array of first-level penalized coefficients indexed by penalty values}
#' \item{gammas}{3-dimensional array of first-level non-penalized coefficients indexed by penalty values}
#' \item{alpha0}{matrix of second-level intercepts indexed by penalty values}
#' \item{alphas}{3-dimensional array of second-level external data coefficients indexed by penalty values}
#' \item{penalty}{vector of first-level penalty values}
#' \item{penalty_ext}{vector of second-level penalty values}
#' \item{penalty_type}{type of penalty applied to first-level predictors}
#' \item{quantile}{quantile for penalty on first-level predictors}
#' \item{penalty_type_ext}{type of penalty applied to second-level external data}
#' \item{quantile_ext}{quantile for penalty on second-level external data}
#' \item{penalty_ratio}{ratio between minimum and maximum penalty (predictors)}
#' \item{penalty_ratio_ext}{ratio between minimum and maximum penalty (external data)}
#' \item{deviance}{fraction of deviance explaned, \eqn{2 x (loglike(saturated model) - loglike(model))}}
#' \item{nlp}{total number of passes over data}
#' \item{custom_mult}{vector of variable-specific penalty multipliers for predictors}
#' \item{custom_mult_ext}{vector of variable-specific penalty multipliers for external data}
#' @importFrom bigmemory is.big.matrix

#' @export
hierr <- function(x,
                  y,
                  external = NULL,
                  unpen = NULL,
                  family = c("gaussian", "binomial"),
                  penalty = definePenalty(),
                  weights = NULL,
                  standardize = c(TRUE, TRUE),
                  intercept = c(TRUE, FALSE),
                  control = list()) {

    # function call
    this.call <- match.call()

    # check error distribution for y
    family <- match.arg(family)

    ## Prepare x and y ##

    # check type of x matrix
    if (is.big.matrix(x)) {
        if (!(bigmemory::describe(x)@description$type %in% c("integer", "double")))
            stop("Error: x contains non-numeric values")
    } else if ("dgCMatrix" %in% class(x)) {
        if (!(typeof(x@x) %in% c("integer", "double")))
            stop("Error: x contains non-numeric values")
    } else {
        stop("Error: x must be a big.matrix, filebacked.big.matrix, or dgCMatrix")
    }

    # check type of y
    y <- as.double(y)

    # check dimensions of x and y
    nr_x <- NROW(x)
    nc_x <- NCOL(x)

    if (nc_x < 2) {
        stop("Error: x must have at least 2 columns")
    }

    y_len <- NROW(y)

    if (y_len != nr_x) {
        stop(paste("Error: Length of y (", y_len, ") not equal to the number of rows of x (", nr_x, ")", sep = ""))
    }

    ## Prepare external ##
    if (!is.null(external)) {

        # check if external is a sparse matrix
        if (is(x, "sparseMatrix")) {
            is_sparse_ext = TRUE
        } else {
            is_sparse_ext = FALSE
            # convert to matrix
            if (class(external) != "matrix") {
                external <- as.matrix(external)
            }
            if (!(typeof(external) %in% c("double", "integer"))) {
                stop("Error: external contains non-numeric values")
            }
        }

        # check dimensions
        nr_ext <- NROW(external)
        nc_ext <- NCOL(external)

        if (nc_x != nr_ext) {
            stop(paste("Error: Number of columns in x (", nc_x, ") not equal to the number of rows in external (", nr_ext, ")", sep = ""))
        }

    } else {
        external <- matrix(vector("numeric", 0), 0, 0)
        nr_ext <- as.integer(0)
        nc_ext <- as.integer(0)
    }

    ## Prepare unpenalized covariates ##
    if (!is.null(unpen)) {

        # check dimensions
        nc_unpen <- NCOL(unpen)

        if (y_len != NROW(unpen)) {
            stop(paste("Error: Length of y (", y_len, ") not equal to the number of rows of unpen (", NROW(unpen), ")", sep = ""))
        }

        # convert unpen to matrix
        if (class(unpen) != "matrix") {
            unpen <- as.matrix(unpen)
        }
        if (!(typeof(x) %in% c("double", "integer"))) {
            stop("Error: unpen contains non-numeric values")
        }
    } else {
        unpen <- matrix(vector("numeric", 0), 0, 0)
        nc_unpen <- as.integer(0)
    }

    # set weights
    if (is.null(weights)) {
        weights <- as.double(rep(1, nr_x))
    } else if (length(weights) != y_len) {
        stop(paste("Error: Length of weights (", length(weights), ") not equal to length of y (", y_len, ")", sep = ""))
    } else if (any(weights < 0)) {
        stop("Error: weights can only contain non-negative values")
    } else {
        weights <- as.double(weights)
    }

    # check penalty object
    penalty <- initialize_penalty(penalty, nr_x, nc_x, nc_unpen, nr_ext, nc_ext, intercept)

    # check control object
    control <- do.call("hierr.control", control)
    control <- initialize_control(control, nc_x, nc_unpen, nc_ext, intercept)

    # fit model
    fit <- fit_model(x = x,
                     y = y,
                     external = external,
                     fixed = unpen,
                     weights_user = weights,
                     intercept = intercept,
                     standardize = standardize,
                     penalty_type = penalty$ptype,
                     cmult = penalty$cmult,
                     quantiles = c(penalty$tau, penalty$tau_ext),
                     num_penalty = c(penalty$num_penalty, penalty$num_penalty_ext),
                     penalty_ratio = c(penalty$penalty_ratio, penalty$penalty_ratio_ext),
                     user_penalty = penalty$user_penalty,
                     user_penalty_ext = penalty$user_penalty_ext,
                     lower_cl = control$lower_limits,
                     upper_cl = control$upper_limits,
                     family = family,
                     thresh = control$tolerance,
                     maxit = control$max_iterations,
                     dfmax = control$dfmax,
                     pmax = control$pmax)

    # check status of model fit
    if (fit$status == 0) {
        # Create arrays ordering coefficients by 1st level penalty / 2nd level penalty
        fit$beta0 <- matrix(fit$beta0, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
        fit$betas <- aperm(array(t(fit$betas), c(penalty$num_penalty_ext, penalty$num_penalty, nc_x)), c(3, 2, 1))
        #fit$deviance <- matrix(fit$deviance, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
        #fit$custom_mult <- penalty$custom_multiplier

        if (intercept[2]) {
            fit$alpha0 <- matrix(fit$alpha0, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
        } else {
            fit$alpha0 <- NULL
        }

        if (nc_ext > 0) {
            #fit$custom_mult_ext <- penalty$custom_multiplier_ext
            fit$alphas <- aperm(array(t(fit$alphas), c(penalty$num_penalty_ext, penalty$num_penalty, nc_ext)), c(3, 2, 1))
            #fit$nzero_alphas <- matrix(fit$nzero_alphas, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
        } else {
            fit$alphas <- NULL
            #fit$nzero_alphas <- NULL
            #fit$penalty_type_ext <- NULL
            #fit$quantile_ext <- NULL
            #fit$penalty_ext <- NULL
            #fit$penalty_ratio_ext <- NULL
        }

        if (nc_unpen > 0) {
            fit$gammas <- aperm(array(t(fit$gammas), c(penalty$num_penalty_ext, penalty$num_penalty, nc_unpen)), c(3, 2, 1))
        } else {
            fit$gammas <- NULL
        }

        #fit$nzero_betas <- matrix(fit$nzero_betas, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
        #fit$num_passes <- matrix(fit$num_passes, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
    } else {
        if (fit$status == -10000) {
            fit$errmsg <- "max iterations reached"
        }
    }

    fit$call <- this.call
    class(fit) <- "hierr"
    return(fit)
}

#' Control function for hierr fitting
#'
#' @description Control function for \code{\link{hierr}} fitting.
#'
#' @param tolerance positive convergence criterion. Default is 1e-07.
#' @param max_iterations maximum number of iterations to run coordinate gradient descent across all penalties before returning an error. Default is 1e+05.
#' @param earlyStop indicator for whether stopping criterion on penalty path based on deviance (i.e. no change in deviance). Default is FALSE.
#' @param dfmax maximum number of variables allowed in model. Default is \eqn{ncol(x) + ncol(external) + intercept[1] + intercept[2]}.
#' @param pmax maximum number of variables with nonzero coefficient estimate. Default is \eqn{min(2 * dfmax + 20, ncol(x) + ncol(external) + intercept[2])}.
#' @param lower_limits vector of lower limits for each coefficient. Default is -Inf for all variables.
#' @param upper_limits vector of upper limits for each coefficient. Default is Inf for all variables.

#' @export
hierr.control <- function(tolerance = 1e-07,
                          max_iterations = 1e+05,
                          earlyStop = FALSE,
                          dfmax = NULL,
                          pmax = NULL,
                          lower_limits = NULL,
                          upper_limits = NULL) {

    if (tolerance < 0) {
        stop("Error: tolerance must be greater than 0")
    }

    if (max_iterations < 0) {
        stop("Error: max_iterations must be a postive integer")
    }

    control_obj <- list(tolerance = tolerance,
                        max_iterations = max_iterations,
                        earlyStop = earlyStop,
                        dfmax = dfmax,
                        pmax = pmax,
                        lower_limits = lower_limits,
                        upper_limits = upper_limits)
}


initialize_penalty <- function(penalty_obj,
                               nr_x,
                               nc_x,
                               nc_unpen,
                               nr_ext,
                               nc_ext,
                               intercept) {
    # check penalty object for x
    if (length(penalty_obj$penalty_type) > 1) {
        if (length(penalty_obj$penalty_type) != nc_x) {
            stop("Error: Length of penalty_type (", length(penalty_obj$penalty_type),") not equal to number of columns in x (", nc_x, ")")
        }
    } else {
        penalty_obj$penalty_type <- rep(penalty_obj$penalty_type, nc_x)
    }

    if (penalty_obj$user_penalty == 0 && is.null(penalty_obj$penalty_ratio)) {
        if (nr_x > nc_x) {
            penalty_obj$penalty_ratio <- 1e-04
        } else {
            penalty_obj$penalty_ratio <- 0.01
        }

        if (penalty_obj$num_penalty < 3) {
            penalty_obj$num_penalty <- 3
            stop("Warning: num_penalty must be at least 3 when automatically computing penalty path")
        }
    }

    if (is.null(penalty_obj$custom_multiplier)) {
        penalty_obj$custom_multiplier <- rep(1.0, nc_x)
    } else if (length(penalty_obj$custom_multiplier) != nc_x) {
        stop("Error: Length of custom_multiplier (", length(penalty_obj$custom_multiplier),") not equal to number of columns in x (", nc_x, ")")
    }

    # check penalty object for external
    if (nc_ext > 0) {
        if (length(penalty_obj$penalty_type_ext) > 1) {
            if (length(penalty_obj$penalty_type_ext) != nc_ext) {
                stop("Error: Length of penalty_type_ext (", length(penalty_obj$penalty_type_ext),") not equal to number of columns in external (", nc_ext, ")")
            }
        } else {
            penalty_obj$penalty_type_ext <- rep(penalty_obj$penalty_type_ext, nc_ext)
        }

        if (penalty_obj$user_penalty_ext == 0 && is.null(penalty_obj$penalty_ratio_ext)) {
            if (nr_ext > nc_ext) {
                penalty_obj$penalty_ratio_ext <- 1e-04
            } else {
                penalty_obj$penalty_ratio_ext <- 0.01
            }

            if (penalty_obj$num_penalty_ext < 3) {
                penalty_obj$num_penalty_ext <- 3
                stop("Warning: num_penalty_ext must be at least 3 when automatically computing penalty path")
            }
        }

        if (is.null(penalty_obj$custom_multiplier_ext)) {
            penalty_obj$custom_multiplier_ext <- rep(1.0, nc_ext)
        } else if (length(penalty_obj$custom_multiplier_ext) != nc_ext && !is.null(external)) {
            stop("Error: Length of custom_multiplier_ext (", length(penalty_obj$custom_multiplier_ext),") not equal to number of columns in external (", nc_ext, ")")
        }
    } else {
        penalty_obj$penalty_type_ext <- NULL
        penalty_obj$num_penalty_ext <- 1
        penalty_obj$penalty_ratio_ext <- 0
        penalty_obj$custom_multiplier_ext <- numeric(0)
    }

    # vectors holding penalty type and multipliers across all variables
    if (intercept[2]) {
        penalty_obj$ptype <- c(penalty_obj$penalty_type, rep(0.0, nc_unpen), 0.0, penalty_obj$penalty_type_ext)
        penalty_obj$cmult <- c(penalty_obj$custom_multiplier, rep(0.0, nc_unpen), 0.0, penalty_obj$custom_multiplier_ext)
    } else {
        penalty_obj$ptype <- c(penalty_obj$penalty_type, rep(0.0, nc_unpen), penalty_obj$penalty_type_ext)
        penalty_obj$cmult <- c(penalty_obj$custom_multiplier, rep(0.0, nc_unpen), penalty_obj$custom_multiplier_ext)
    }
    return(penalty_obj)
}

initialize_control <- function(control_obj, nc_x, nc_unpen, nc_ext, intercept) {
    if (is.null(control_obj$dfmax)) {
        control_obj$dfmax <- as.integer(nc_x + nc_ext + nc_unpen + intercept[1] + intercept[2])
    } else if (control_obj$dfmax < 0) {
        stop("Error: dfmax can only contain postive integers")
    }

    if (is.null(control_obj$pmax)) {
        control_obj$pmax <- as.integer(min(2 * control_obj$dfmax + 20, nc_x + nc_ext + nc_unpen + intercept[2]))
    } else if (control_obj$pmax < 0) {
        stop("Error: pmax can only contain positive integers")
    }

    if (is.null(control_obj$lower_limits)) {
        control_obj$lower_limits <- rep(-Inf, nc_x + nc_ext + nc_unpen + intercept[2])
    } else if (length(control_obj$lower_limits) != nc_x + nc_ext + nc_unpen) {
        stop("Error: Length of lower_limits (", length(control_obj$lower_limits), ") not equal to sum of number of columns in x, unpen, and external (", nc_x + nc_ext + nc_unpen, ")")
    } else if (intercept[2]) {
        control_obj$lower_limits <- c(control_obj$lower_limits[1:(nc_x + nc_unpen)],
                                      -Inf,
                                      control_obj$lower_limits[(nc_x + nc_unpen + 1):length(control_obj$lower_limits)])
    }

    if (is.null(control_obj$upper_limits)) {
        control_obj$upper_limits <- rep(Inf, nc_x + nc_ext + nc_unpen + intercept[2])
    } else if (length(control_obj$upper_limits) != nc_x + nc_ext + nc_unpen) {
        stop("Error: Length of upper_limits (", length(control_obj$upper_limits), ") not equal to sum of number of columns in x, unpen, and external (", nc_x + nc_ext + nc_unpen, ")")
    } else if (intercept[2]) {
        control_obj$upper_limits <- c(control_obj$upper_limits[1:(nc_x + nc_unpen)],
                                      -Inf,
                                      control_obj$upper_limits[(nc_x + nc_unpen + 1):length(control_obj$upper_limits)])
    }
    return(control_obj)
}
