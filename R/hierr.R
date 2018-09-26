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
#' @param x predictor design matrix of dimension \eqn{n x p}
#' @param y outcome vector of length \eqn{n}
#' @param external (optional) external data design matrix of dimension \eqn{p x q}. Can be class "matrix" or "dgCMatrix".
#' @param unpen (optional) unpenalized predictor design matrix
#' @param family error distribution for outcome variable
#' @param penalty specifies regularization object for x and external. See \code{\link{definePenalty}} for more details.
#' @param weights optional vector of observation-specific weights. Default is 1 for all observations.
#' @param standardize indicates whether x and/or external should be standardized. Default is c(TRUE, TRUE).
#' @param intercept indicates whether an intercept term is included for x and/or external. Default is c(TRUE, TRUE).
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


#' @export
hierr <- function(x,
                  y,
                  external = NULL,
                  unpen = NULL,
                  family = c("gaussian"),
                  penalty = definePenalty(),
                  weights = NULL,
                  standardize = c(TRUE, TRUE),
                  intercept = c(TRUE, TRUE),
                  control = list()) {

    # function call
    this.call <- match.call()

    # check error distribution for y
    family <- match.arg(family)

    ## Prepare x and y ##

    # check if x is a sparse matrix
    if (inherits(x, "sparseMatrix")) {
        stop("Error: sparse matrices only supported for external data")
    } else {
        # convert x to matrix
        if (class(x) != "matrix") {
            x <- as.matrix(x)
        }
        if (!(typeof(x) %in% c("double", "integer"))) {
            stop("Error: x contains non-numeric values")
        }
    }

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
    is_sparse = FALSE
    if (!is.null(external)) {

        # check if external is a sparse matrix
        if (class(external) %in% "dgCMatrix") {
            is_sparse = TRUE
        } else {
            # convert to matrix
            if (class(external) != "matrix") {
                external <- as.matrix(external)
            }
            if (!(typeof(x) %in% c("double", "integer"))) {
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
        external <- matrix(numeric(0), nrow = 0, ncol = 0)
        nr_ext <- as.integer(0)
        nc_ext <- as.integer(0)
    }

    ## Prepare unpenalized covariates ##
    if (!is.null(unpen)) {

        # check dimensions
        nr_unpen <- NROW(unpen)
        nc_unpen <- NCOL(unpen)

        if (y_len != nr_unpen) {
            stop(paste("Error: Length of y (", y_len, ") not equal to the number of rows of unpen (", nr_unpen, ")", sep = ""))
        }

        # convert unpen to matrix
        if (class(unpen) != "matrix") {
            unpen <- as.matrix(unpen)
        }
        if (!(typeof(x) %in% c("double", "integer"))) {
            stop("Error: unpen contains non-numeric values")
        }
    } else {
        unpen <- matrix(numeric(0), nrow = 0, ncol = 0)
        nr_unpen <- as.integer(0)
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

    # check penalty object for x
    if (length(penalty$penalty_type) > 1) {
        if (length(penalty$penalty_type) != nc_x) {
            stop("Error: Length of penalty_type (", length(penalty$penalty_type),") not equal to number of columns in x (", nc_x, ")")
        }
    } else {
        penalty$penalty_type <- rep(penalty$penalty_type, nc_x)
    }

    if (penalty$user_penalty == 0 && is.null(penalty$penalty_ratio)) {
        if (nr_x > nc_x) {
            penalty$penalty_ratio <- 1e-04
        } else {
            penalty$penalty_ratio <- 0.01
        }

        if (penalty$num_penalty < 3) {
            penalty$num_penalty <- 3
            stop("Warning: num_penalty must be at least 3 when automatically computing penalty path")
        }
    }

    if (is.null(penalty$custom_multiplier)) {
        penalty$custom_multiplier <- rep(1.0, nc_x)
    } else if (length(penalty$custom_multiplier) != nc_x) {
        stop("Error: Length of custom_multiplier (", length(penalty$custom_multiplier),") not equal to number of columns in x (", nc_x, ")")
    }

    # check penalty object for external
    if (nc_ext > 0) {
        if (length(penalty$penalty_type_ext) > 1) {
            if (length(penalty$penalty_type_ext) != nc_ext) {
                stop("Error: Length of penalty_type_ext (", length(penalty$penalty_type_ext),") not equal to number of columns in external (", nc_ext, ")")
            }
        } else {
            penalty$penalty_type_ext <- rep(penalty$penalty_type_ext, nc_ext)
        }

        if (penalty$user_penalty_ext == 0 && is.null(penalty$penalty_ratio_ext)) {
            if (nr_ext > nc_ext) {
                penalty$penalty_ratio_ext <- 1e-04
            } else {
                penalty$penalty_ratio_ext <- 0.01
            }

            if (penalty$num_penalty_ext < 3) {
                penalty$num_penalty_ext <- 3
                stop("Warning: num_penalty_ext must be at least 3 when automatically computing penalty path")
            }
        }

        if (is.null(penalty$custom_multiplier_ext)) {
            penalty$custom_multiplier_ext <- rep(1.0, nc_ext)
        } else if (length(penalty$custom_multiplier_ext) != nc_ext && !is.null(external)) {
            stop("Error: Length of custom_multiplier_ext (", length(penalty$custom_multiplier_ext),") not equal to number of columns in external (", nc_ext, ")")
        }
    } else {
        penalty$penalty_type_ext <- NULL
        penalty$num_penalty_ext <- 1
        penalty$penalty_ratio_ext <- 0
        penalty$custom_multiplier_ext <- numeric(0)
    }

    # vectors holding penalty type and multipliers across all variables
    if (intercept[2]) {
        penalty$ptype <- c(penalty$penalty_type, rep(0.0, nc_unpen), 0.0, penalty$penalty_type_ext)
        penalty$cmult <- c(penalty$custom_multiplier, rep(0.0, nc_unpen), 0.0, penalty$custom_multiplier_ext)
    } else {
        penalty$ptype <- c(penalty$penalty_type, rep(0.0, nc_unpen), penalty$penalty_type_ext)
        penalty$cmult <- c(penalty$custom_multiplier, rep(0.0, nc_unpen), penalty$custom_multiplier_ext)
    }

    # check control object
    control <- do.call("hierr.control", control)

    if (is.null(control$dfmax)) {
        control$dfmax <- as.integer(nc_x + nc_ext + nc_unpen + intercept[1] + intercept[2])
    } else if (control$dfmax < 0) {
        stop("Error: dfmax can only contain postive integers")
    }

    if (is.null(control$pmax)) {
        control$pmax <- as.integer(min(2 * control$dfmax + 20, nc_x + nc_ext + nc_unpen + intercept[2]))
    } else if (control$pmax < 0) {
        stop("Error: pmax can only contain positive integers")
    }

    if (is.null(control$lower_limits)) {
        control$lower_limits <- rep(-Inf, nc_x + nc_ext + nc_unpen + intercept[2])
    } else if (length(control$lower_limits) != nc_x + nc_ext + nc_unpen) {
        stop("Error: Length of lower_limits (", length(control$lower_limits), ") not equal to sum of number of columns in x, unpen, and external (", nc_x + nc_ext + nc_unpen, ")")
    } else if (intercept[2]) {
        control$lower_limits <- c(control$lower_limits[1:(nc_x + nc_unpen)], -Inf, control$lower_limits[(nc_x + nc_unpen + 1):length(control$lower_limits)])
    }

    if (is.null(control$upper_limits)) {
        control$upper_limits <- rep(Inf, nc_x + nc_ext + nc_unpen + intercept[2])
    } else if (length(control$upper_limits) != nc_x + nc_ext + nc_unpen) {
        stop("Error: Length of upper_limits (", length(control$upper_limits), ") not equal to sum of number of columns in x, unpen, and external (", nc_x + nc_ext + nc_unpen, ")")
    } else if (intercept[2]) {
        control$upper_limits <- c(control$upper_limits[1:(nc_x + nc_unpen)], -Inf, control$upper_limits[(nc_x + nc_unpen + 1):length(control$upper_limits)])
    }

    # fit model based on distribution
    fit <- do.call(family, list(x = x,
                                y = y,
                                external = external,
                                unpen = unpen,
                                weights = weights,
                                penalty = penalty,
                                isd = standardize,
                                intr = intercept,
                                is_sparse = is_sparse,
                                control = control))


    # check status of model fit
    if (fit$status == 0) {
        # Create arrays ordering coefficients by 1st level penalty / 2nd level penalty
        fit$beta0 <- matrix(fit$beta0, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
        fit$betas <- aperm(array(t(fit$betas), c(penalty$num_penalty_ext, penalty$num_penalty, nc_x)), c(3, 2, 1))
        fit$deviance <- matrix(fit$deviance, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
        fit$custom_mult <- penalty$custom_multiplier

        if (intercept[2]) {
            fit$alpha0 <- matrix(fit$alpha0, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
        } else {
            fit$alpha0 <- NULL
        }

        if (nc_ext > 0) {
            fit$custom_mult_ext <- penalty$custom_multiplier_ext
            fit$alphas <- aperm(array(t(fit$alphas), c(penalty$num_penalty_ext, penalty$num_penalty, nc_ext)), c(3, 2, 1))
            fit$nzero_alphas <- matrix(fit$nzero_alphas, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
        } else {
            fit$alphas <- NULL
            fit$nzero_alphas <- NULL
            fit$penalty_type_ext <- NULL
            fit$quantile_ext <- NULL
            fit$penalty_ext <- NULL
            fit$penalty_ratio_ext <- NULL
        }

        if (nc_unpen > 0) {
            fit$gammas <- aperm(array(t(fit$gammas), c(penalty$num_penalty_ext, penalty$num_penalty, nc_unpen)), c(3, 2, 1))
        } else {
            fit$gammas <- NULL
        }

        fit$nzero_betas <- matrix(fit$nzero_betas, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
        fit$num_passes <- matrix(fit$num_passes, nrow = penalty$num_penalty, ncol = penalty$num_penalty_ext, byrow = TRUE)
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
