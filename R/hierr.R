#' @useDynLib hierr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Fit hierarchical regularized regression model
#'
#' @param x predictor design matrix of dimension n x p
#' @param y outcome vector of length n
#' @param external external data design matrix of dimension p x q
#' @param family error distribution for outcome variable
#' @param penalty specifies regularization object for x and external. See \code{\link{definePenalty}} for more details.
#' @param weights optional vector of observation-specific weights. Default is 1 for all observations.
#' @param standardize indicates whether x and/or external should be standardized. Default is c(TRUE, TRUE).
#' @param intercept indicates whether an intercept term is included for x and/or external. Default is c(TRUE, TRUE).
#' @param control specifies hierr control object. See \code{\link{hierr.control}} for more details.

#' @export
hierr <- function(x,
                  y,
                  external,
                  family = c("gaussian", "binomial", "poisson"),
                  penalty = definePenalty(0, 1),
                  weights = NULL,
                  standardize = c(TRUE, TRUE),
                  intercept = c(TRUE, TRUE),
                  control = list()) {

    # function call
    this.call <- match.call()

    # check error distribution for y
    family <- match.arg(family)

    nr_x <- nrow(x)
    nc_x <- ncol(x)
    nr_ext <- nrow(external)
    nc_ext <- ncol(external)

    # check penalty object
    if (penalty$user_penalty == 0 && is.null(penalty$penalty_ratio)) {
        if (nr_x > nc_x) {
            penalty$penalty_ratio <- 1e-04
        } else {
            penalty$penalty_ratio <- 0.01
        }
    }

    if (penalty$user_penalty_ext == 0 && is.null(penalty$penalty_ratio_ext)) {
        if (nr_ext > nc_ext) {
            penalty$penalty_ratio_ext <- 1e-04
        } else {
            penalty$penalty_ratio_ext <- 0.01
        }
    }

    if (is.null(penalty$custom_multiplier)) {
        penalty$custom_multiplier <- as.double(rep(1, nc_x))
    } else {
        if (length(penalty$custom_multiplier) != nc_x) {
            stop("Length of custom_multiplier (", length(penalty$custom_multiplier),") not equal to number of columns of x (", nc_x, ")")
        }
    }

    if (is.null(penalty$custom_multiplier_ext)) {
        penalty$custom_multiplier_ext <- as.double(rep(1, nc_ext))
    } else {
        if (length(penalty$custom_multiplier_ext) != nc_ext) {
            stop("Length of custom_multiplier_ext (", length(penalty$custom_multiplier_ext),") not equal to number of columns of external (", nc_ext, ")")
        }
    }

    # check dimensions
    if (nc_x < 2 || nc_ext < 2) {
        stop("Both x and external must have at least 2 columns")
    }

    y_len <- ifelse(is.null(dim(y)), length(y), dim(y)[1])

    if (y_len != nr_x) {
        stop(paste("Number of observations in y (", y_len, ") not equal to the number of rows of x (", nr_x, ")", sep = ""))
    }

    if (nc_x != nr_ext) {
        stop(paste("Number of columns in x (", nc_x, ") not equal to the number of rows in external (", nr_ext, ")", sep = ""))
    }

    # set weights
    if (is.null(weights)) {
        weights <- as.double(rep(1, nr_x))
    } else if (length(weights) != nr_x) {
        stop(paste("Number of elements in weights (", length(weights), ") not equal to the number of rows of x (", nr_x, ")", sep = ""))
    } else if (any(weights) < 0) {
        stop("weights can only contain non-negative values")
    } else {
        weights <- as.double(weights)
    }

    # convert to matrices
    if (!(class(x)) != "matrix") {
        x <- as.matrix(x)
    }
    if (typeof(x) != "double") {
        stop("x contains non-numeric values")
    }

    if (!(class(external)) != "matrix") {
        external <- as.matrix(external)
    }
    if (typeof(external) != "double") {
        stop("external contains non-numeric values")
    }

    # check control object
    control <- do.call("hierr.control", control)

    if (is.null(control$dfmax)) {
        control$dfmax <- as.integer(nc_x + nc_ext + intercept[1] + intercept[2])
    } else {
        control$dfmax <- as.integer(control$dfmax)
    }

    if (is.null(control$pmax)) {
        control$pmax <- as.integer(min(2 * control$dfmax + 20, nc_x + nc_ext + 1))
    } else {
        control$pmax <- as.integer(control$pmax)
    }

    if (is.null(control$lower_limits)) {
        control$lower_limits <- rep(-Inf, nc_x + nc_ext)
    } else if (length(control$lower_limits) != nc_x + nc_ext) {
        stop("Length of lower_limits (", length(control$lower_limits), ") not equal to sum of number of columns in x and external (", nc_x + nc_ext, ")")
    }

    if (is.null(control$upper_limits)) {
        control$upper_limits <- rep(Inf, nc_x + nc_ext)
    } else if (length(control$upper_limits) != nc_x + nc_ext) {
        stop("Length of upper_limits (", length(control$upper_limits), ") not equal to sum of number of columns in x and external (", nc_x + nc_ext, ")")
    }

    fit <- do.call(family, list(x = x,
                                y = y,
                                external = external,
                                weights = weights,
                                penalty = penalty,
                                isd = standardize,
                                intr = intercept,
                                control = control)
                   )
}

#' Control function for hierr fitting
#'
#'@description Control function for \code{\link{hierr}} fitting.
#'
#' @param tolerance positive convergence criterion. Default is 1e-08.
#' @param max_iterations maximum number of iterations to run coordinate gradient descent across all penalties before returning an error. Default is 1e+05.
#' @param dfmax maximum number of variables allowed in model. Default is \eqn{p + 1}.
#' @param pmax maximum number of variables with nonzero coefficient estimate. Default is \eqn{min(2*dfmax + 20, p)}.
#' @param lower_limits vector of lower limits for each coefficient. Default is -Inf.
#' @param upper_limits vector of upper limits for each coefficient. Default is Inf.

hierr.control <- function(tolerance = 1e-08,
                          max_iterations = 1e+05,
                          dfmax = NULL,
                          pmax = NULL,
                          lower_limits = NULL,
                          upper_limits = NULL) {

    if (tolerance < 0) {
        stop("tolerance must be greater than 0")
    }

    if (max_iterations < 0) {
        stop("max_iterations must be greater than 0")
    }

    structure(list(tolerance = tolerance,
                   max_iterations = max_iterations,
                   dfmax = dfmax,
                   pmax = pmax,
                   lower_limits = lower_limits,
                   upper_limits = upper_limits)
              )
}
