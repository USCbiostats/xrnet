#' @useDynLib hierr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Fit hierarchical regularized regression model
#'
#' @param x predictor design matrix of dimension \eqn{n x p}
#' @param y outcome vector of length \eqn{n}
#' @param external (optional) external data design matrix of dimension \eqn{p x q}
#' @param family error distribution for outcome variable
#' @param penalty specifies regularization object for x and external. See \code{\link{definePenalty}} for more details.
#' @param weights optional vector of observation-specific weights. Default is 1 for all observations.
#' @param standardize indicates whether x and/or external should be standardized. Default is c(TRUE, TRUE).
#' @param intercept indicates whether an intercept term is included for x and/or external. Default is c(TRUE, TRUE).
#' @param control specifies hierr control object. See \code{\link{hierr.control}} for more details.

#' @export
hierr <- function(x,
                  y,
                  external = NULL,
                  family = c("gaussian"),
                  penalty = definePenalty(0, 1),
                  weights = NULL,
                  standardize = c(TRUE, TRUE),
                  intercept = c(TRUE, TRUE),
                  control = list()) {

    # function call
    this.call <- match.call()

    # check error distribution for y
    family <- match.arg(family)

    ## Prepare x and y ##

    # check dimensions of x and y
    nr_x <- nrow(x)
    nc_x <- ncol(x)

    if (nc_x < 2) {
        stop("Error: x must have at least 2 columns")
    }

    y_len <- ifelse(is.null(dim(y)), length(y), dim(y)[1])

    if (y_len != nr_x) {
        stop(paste("Error: Length of y (", y_len, ") not equal to the number of rows of x (", nr_x, ")", sep = ""))
    }

    # set weights
    if (is.null(weights)) {
        weights <- as.double(rep(1, nr_x))
    } else if (length(weights) != nr_x) {
        stop(paste("Error: Length of weights (", length(weights), ") not equal to the number of rows of x (", nr_x, ")", sep = ""))
    } else if (any(weights < 0)) {
        stop("Error: weights can only contain non-negative values")
    } else {
        weights <- as.double(weights)
    }

    # convert x to matrix
    if (!(class(x)) != "matrix") {
        x <- as.matrix(x)
    }
    if (typeof(x) != "double") {
        stop("Error: x contains non-numeric values")
    }

    ## Prepare external ##
    ext_exist = FALSE
    if (!is.null(external)) {
        ext_exist = TRUE

        # check dimensions
        nr_ext <- nrow(external)
        nc_ext <- ncol(external)

        if (nc_x != nr_ext) {
            stop(paste("Error: Number of columns in x (", nc_x, ") not equal to the number of rows in external (", nr_ext, ")", sep = ""))
        }

        # convert to matrix
        if (!(class(external)) != "matrix") {
            external <- as.matrix(external)
        }
        if (typeof(external) != "double") {
            stop("Error: external contains non-numeric values")
        }
    } else {
        external <- matrix(numeric(0), nrow = 0, ncol = 0)
        nr_ext <- as.integer(0)
        nc_ext <- as.integer(0)
    }

    # check penalty object for x
    if (penalty$user_penalty == 0 && is.null(penalty$penalty_ratio)) {
        if (nr_x > nc_x) {
            penalty$penalty_ratio <- 1e-04
        } else {
            penalty$penalty_ratio <- 0.01
        }
    }

    if (is.null(penalty$custom_multiplier)) {
        penalty$custom_multiplier <- rep(1.0, nc_x)
    } else if (length(penalty$custom_multiplier) != nc_x) {
        stop("Error: Length of custom_multiplier (", length(penalty$custom_multiplier),") not equal to number of columns of x (", nc_x, ")")
    }

    # check penalty object for ext
    if (ext_exist) {
        if (penalty$user_penalty_ext == 0 && is.null(penalty$penalty_ratio_ext)) {
            if (nr_ext > nc_ext) {
                penalty$penalty_ratio_ext <- 1e-04
            } else {
                penalty$penalty_ratio_ext <- 0.01
            }
        }

        if (is.null(penalty$custom_multiplier_ext)) {
            penalty$custom_multiplier_ext <- rep(1.0, nc_ext)
        } else if (length(penalty$custom_multiplier_ext) != nc_ext && !is.null(external)) {
            stop("Error: Length of custom_multiplier_ext (", length(penalty$custom_multiplier_ext),") not equal to number of columns of external (", nc_ext, ")")
        }
    } else {
        penalty$num_penalty_ext <- 1
        penalty$penalty_ratio_ext <- 0
        penalty$custom_multiplier_ext <- numeric(0)
    }

    if (intercept[2]) {
        penalty$cmult <- c(penalty$custom_multiplier, 0.0, penalty$custom_multiplier_ext)
    } else {
        penalty$cmult <- c(penalty$custom_multiplier, penalty$custom_multiplier_ext)
    }
    penalty$ptype <- c(rep(penalty$penalty_type, nc_x), rep(penalty$penalty_type_ext, nc_ext + intercept[2]))

    # check control object
    control <- do.call("hierr.control", control)

    if (is.null(control$dfmax)) {
        control$dfmax <- as.integer(nc_x + nc_ext + intercept[1] + intercept[2])
    } else if (control$dfmax < 0) {
        stop("Error: dfmax can only contain postive integers")
    }

    if (is.null(control$pmax)) {
        control$pmax <- as.integer(min(2 * control$dfmax + 20, nc_x + nc_ext + intercept[2]))
    } else if (control$pmax < 0) {
        stop("Error: pmax can only contain positive integers")
    }

    if (is.null(control$lower_limits)) {
        control$lower_limits <- rep(-Inf, nc_x + nc_ext + intercept[2])
    } else if (length(control$lower_limits) != nc_x + nc_ext) {
        stop("Error: Length of lower_limits (", length(control$lower_limits), ") not equal to sum of number of columns in x and external (", nc_x + nc_ext, ")")
    } else if (intercept[2]) {
        control$lower_limits <- c(control$lower_limits[1:nc_x], -Inf, control$lower_limits[(nc_x + 1):(nc_x + nc_ext)])
    }

    if (is.null(control$upper_limits)) {
        control$upper_limits <- rep(Inf, nc_x + nc_ext + intercept[2])
    } else if (length(control$upper_limits) != nc_x + nc_ext) {
        stop("Error: Length of upper_limits (", length(control$upper_limits), ") not equal to sum of number of columns in x and external (", nc_x + nc_ext, ")")
    } else if (intercept[2]) {
        control$upper_limits <- c(control$upper_limits[1:nc_x], -Inf, control$upper_limits[(nc_x + 1):(nc_x + nc_ext)])
    }

    # fit model based on distribution
    fit <- do.call(family, list(x = x,
                                y = y,
                                external = external,
                                weights = weights,
                                penalty = penalty,
                                isd = standardize,
                                intr = intercept,
                                control = control))

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

    if (ext_exist) {
        fit$custom_mult_ext <- penalty$custom_multiplier_ext
        fit$alphas <- aperm(array(t(fit$alphas), c(penalty$num_penalty_ext, penalty$num_penalty, nc_ext)), c(3, 2, 1))
    } else {
        fit$alphas <- NULL
        fit$penalty_type_ext <- NULL
        fit$penalty_ext <- NULL
        fit$penalty_ratio_ext <- NULL
    }

    fit$call <- this.call
    class(fit) <- "hierr"
    return(fit)
}

#' Control function for hierr fitting
#'
#' @description Control function for \code{\link{hierr}} fitting.
#'
#' @param tolerance positive convergence criterion. Default is 1e-08.
#' @param max_iterations maximum number of iterations to run coordinate gradient descent across all penalties before returning an error. Default is 1e+05.
#' @param dfmax maximum number of variables allowed in model. Default is \eqn{ncol(x) + ncol(external) + intercept[1] + intercept[2]}.
#' @param pmax maximum number of variables with nonzero coefficient estimate. Default is \eqn{min(2 * dfmax + 20, ncol(x) + ncol(external) + intercept[2])}.
#' @param lower_limits vector of lower limits for each coefficient. Default is -Inf for all variables.
#' @param upper_limits vector of upper limits for each coefficient. Default is Inf for all variables.

#' @export
hierr.control <- function(tolerance = 1e-08,
                          max_iterations = 1e+05,
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

    structure(list(tolerance = tolerance,
                   max_iterations = max_iterations,
                   dfmax = dfmax,
                   pmax = pmax,
                   lower_limits = lower_limits,
                   upper_limits = upper_limits))
}
