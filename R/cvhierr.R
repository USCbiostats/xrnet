#' k-fold cross-validation for hierarchical regularized regression
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom bigmemory describe
#' @importFrom bigmemory attach.big.matrix
#'
#' @description k-fold cross-validation for hierarchical regularized regression \code{\link{hierr}}
#'
#' @param x predictor design matrix of dimension \eqn{n x p}, matrix options include:
#' \itemize{
#'    \item matrix
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
#' @param loss loss function for cross-validation. Options include:
#' \itemize{
#'    \item mse (Mean Squared Error)
#'    \item mae (Mean Absolute Error)
#'    \item auc (Area under the curve)
#' }
#' @param nfolds number of folds for cross-validation. Default is 5.
#' @param foldid (optional) vector that identifies user-specified fold for each observation. If NULL, folds are automatically generated.
#' @param parallel use \code{foreach} function to fit folds in parallel if TRUE, must register cluster (\code{doParallel}) before using.
#' @param control specifies hierr control object. See \code{\link{hierr.control}} for more details.
#' @return A list of class \code{cvhierr} with components
#' \item{cv_mean}{mean cross-validated error for each penalty combination. Object returned is
#' a vector if there is no external data (external = NULL) and matrix if there is external data.}
#' \item{cv_sd}{estimated standard deviation for cross-validated errors}
#' \item{opt_loss}{the value for the optimal cross-validated error}
#' \item{opt_penalty}{first-level penalty value that achieves the optimal loss}
#' \item{opt_penalty_ext}{second-level penalty value that achieves the optimal loss (if external data is present)}
#' \item{fit_train}{fitted hierr object for all training data, see \code{\link{hierr}} for details of object}

#' @export
cv_hierr <- function(x,
                     y,
                     external = NULL,
                     unpen = NULL,
                     family = c("gaussian", "binomial"),
                     penalty = definePenalty(),
                     weights = NULL,
                     standardize = c(TRUE, TRUE),
                     intercept = c(TRUE, FALSE),
                     loss = c("mse", "mae", "auc"),
                     nfolds = 5,
                     foldid = NULL,
                     parallel = FALSE,
                     control = list())
{

    # Check family argument
    family <- match.arg(family)

    # Set measure used to assess model prediction performance
    if (missing(loss)) {
        if (family == "gaussian")
            loss <- "mse"
        else if (family == "binomial")
            loss <- "auc"
    } else {
        loss <- match.arg(loss)
        loss_available <- TRUE
        if (family == "gaussian" && !(loss %in% c("mse", "mae")))
            loss_available <- FALSE
        else if(family == "binomial" && !(loss %in% c("auc")))
            loss_available <- FALSE
        if (!loss_available)
            stop(paste0("Error: loss = '", loss, "' is not available for family = '", family,"'"))
    }

    # check type of x matrix
    if (is(x, "matrix")) {
        if (!(typeof(x) %in% c("integer", "double")))
            stop("Error: x contains non-numeric values")
        mattype_x <- 1
    }
    else if (is.big.matrix(x)) {
        if (!(bigmemory::describe(x)@description$type %in% c("integer", "double")))
            stop("Error: x contains non-numeric values")
        mattype_x <- 2
    } else if ("dgCMatrix" %in% class(x)) {
        if (!(typeof(x@x) %in% c("integer", "double")))
            stop("Error: x contains non-numeric values")
        mattype_x <- 3
    } else {
        stop("Error: x must be a big.matrix, filebacked.big.matrix, or dgCMatrix")
    }

    # check y type
    y <- drop(as.numeric(y))

    # Get arguments to cvhierr() function and filter for calls to fitting procedure
    hierr_call <- match.call(expand.dots = TRUE)
    cv_args <- match(c("loss", "nfolds", "foldid", "parallel"), names(hierr_call), FALSE)

    if (any(cv_args)) {
        hierr_call <- hierr_call[-cv_args]
    }
    hierr_call[[1]] <- as.name("hierr")

    # Set sample size / weights
    n <- length(y)
    if (is.null(weights)) {
        weights <- rep(1, n)
    }

    # Fit model on all training data
    hierr_object <- hierr(x = x,
                          y = y,
                          external = external,
                          unpen = unpen,
                          family = family,
                          weights = weights,
                          standardize = standardize,
                          intercept = intercept,
                          penalty = penalty,
                          control = control)
    hierr_object$call <- hierr_call

    # Check whether fixed and external are empty
    if (is.null(unpen)) {
        unpen <- matrix(vector("numeric", 0), 0, 0)
        nc_unpen <- as.integer(0)
    } else {
        nc_unpen <- NCOL(unpen)
    }
    if (is.null(external)) {
        external <- matrix(vector("numeric", 0), 0, 0)
        nc_ext <- as.integer(0)
    } else {
        nc_ext <- NCOL(external)
    }

    # Prepare penalty and control object for folds
    penalty_fold <- penalty
    penalty_fold$user_penalty <- hierr_object$penalty
    penalty_fold$user_penalty_ext <- hierr_object$penalty_ext

    penalty_fold <- initialize_penalty(penalty_fold,
                                       NROW(x),
                                       NCOL(x),
                                       nc_unpen,
                                       NROW(external),
                                       nc_ext,
                                       intercept)

    num_pen <- penalty_fold$num_penalty
    num_pen_ext <- penalty_fold$num_penalty_ext

    control <- do.call("hierr.control", control)
    control <- initialize_control(control, NCOL(x), nc_unpen, nc_ext, intercept)

    # Randomly sample observations into folds / check nfolds
    if (is.null(foldid)) {
        if (nfolds < 2)
            stop("number of folds (nfolds) must be at least 2")
        foldid <- sample(rep(seq(nfolds), length = n))
    } else {
        if (length(foldid) != n) {
            stop("Error: length of foldid (", foldid, ") not equal to number of observations (", n, ")")
        }
        foldid <- as.numeric(factor(foldid))
        nfolds <- length(unique(foldid))
        if (nfolds < 2)
            stop("number of folds (nfolds) must be at least 2")
    }

    # Run k-fold CV
    if (parallel) {
        if (is.big.matrix(x)) {
            xdesc <- describe(x)
            errormat <- foreach(k = 1L:nfolds, .packages = c("hierr", "bigmemory"), .combine = cbind) %dopar% {
                weights_train <- weights
                weights_train[foldid == k] <- 0.0
                test_idx <- as.integer(which(foldid == k) - 1)
                xref <- attach.big.matrix(xdesc)

                # Get errors for k-th fold
                error_vec <- fit_model_cv(x = xref,
                                          mattype_x = mattype_x,
                                          y = y,
                                          external = external,
                                          fixed = unpen,
                                          weights_user = weights_train,
                                          intercept = intercept,
                                          standardize = standardize,
                                          penalty_type = penalty_fold$ptype,
                                          cmult = penalty_fold$cmult,
                                          quantiles = c(penalty_fold$quantile, penalty_fold$quantile_ext),
                                          num_penalty = c(penalty_fold$num_penalty, penalty_fold$num_penalty_ext),
                                          penalty_ratio = c(penalty_fold$penalty_ratio, penalty_fold$penalty_ratio_ext),
                                          user_penalty = penalty_fold$user_penalty,
                                          user_penalty_ext = penalty_fold$user_penalty_ext,
                                          lower_cl = control$lower_limits,
                                          upper_cl = control$upper_limits,
                                          family = family,
                                          user_loss = loss,
                                          test_idx = test_idx,
                                          thresh = control$tolerance,
                                          maxit = control$max_iterations,
                                          dfmax = control$dfmax,
                                          pmax = control$pmax)
            }
        } else {
            errormat <- foreach(k = 1L:nfolds, .packages = c("hierr", "Matrix"), .combine = cbind) %dopar% {
                weights_train <- weights
                weights_train[foldid == k] <- 0.0
                test_idx <- as.integer(which(foldid == k) - 1)

                # Get errors for k-th fold
                error_vec <- fit_model_cv(x = x,
                                          mattype_x = mattype_x,
                                          y = y,
                                          external = external,
                                          fixed = unpen,
                                          weights_user = weights_train,
                                          intercept = intercept,
                                          standardize = standardize,
                                          penalty_type = penalty_fold$ptype,
                                          cmult = penalty_fold$cmult,
                                          quantiles = c(penalty_fold$quantile, penalty_fold$quantile_ext),
                                          num_penalty = c(penalty_fold$num_penalty, penalty_fold$num_penalty_ext),
                                          penalty_ratio = c(penalty_fold$penalty_ratio, penalty_fold$penalty_ratio_ext),
                                          user_penalty = penalty_fold$user_penalty,
                                          user_penalty_ext = penalty_fold$user_penalty_ext,
                                          lower_cl = control$lower_limits,
                                          upper_cl = control$upper_limits,
                                          family = family,
                                          user_loss = loss,
                                          test_idx = test_idx,
                                          thresh = control$tolerance,
                                          maxit = control$max_iterations,
                                          dfmax = control$dfmax,
                                          pmax = control$pmax)
            }
        }
    } else {
        errormat <- matrix(NA, nrow = num_pen * num_pen_ext, ncol = nfolds)
        for (k in 1:nfolds) {
            # Split into test and train for k-th fold
            weights_train <- weights
            weights_train[foldid == k] <- 0.0
            test_idx <- as.integer(which(foldid == k) - 1)

            # Fit model on k-th training fold
            errormat[, k] <- fit_model_cv(x = x,
                                          mattype_x = mattype_x,
                                          y = y,
                                          external = external,
                                          fixed = unpen,
                                          weights_user = weights_train,
                                          intercept = intercept,
                                          standardize = standardize,
                                          penalty_type = penalty_fold$ptype,
                                          cmult = penalty_fold$cmult,
                                          quantiles = c(penalty_fold$quantile, penalty_fold$quantile_ext),
                                          num_penalty = c(penalty_fold$num_penalty, penalty_fold$num_penalty_ext),
                                          penalty_ratio = c(penalty_fold$penalty_ratio, penalty_fold$penalty_ratio_ext),
                                          user_penalty = penalty_fold$user_penalty,
                                          user_penalty_ext = penalty_fold$user_penalty_ext,
                                          lower_cl = control$lower_limits,
                                          upper_cl = control$upper_limits,
                                          family = family,
                                          user_loss = loss,
                                          test_idx = test_idx,
                                          thresh = control$tolerance,
                                          maxit = control$max_iterations,
                                          dfmax = control$dfmax,
                                          pmax = control$pmax)
        }
    }
    cv_mean <- rowMeans(errormat)
    cv_sd <- sqrt(rowSums((errormat - cv_mean)^2) / (nfolds - 1))
    cv_mean <- matrix(cv_mean, nrow = num_pen, byrow = TRUE)
    cv_sd <- matrix(cv_sd, nrow = num_pen, byrow = TRUE)
    row.names(cv_mean) <- rev(sort(hierr_object$penalty))
    row.names(cv_sd) <- rev(sort(hierr_object$penalty))
    if (num_pen_ext > 1) {
        colnames(cv_mean) <- rev(sort(hierr_object$penalty_ext))
        colnames(cv_sd) <- rev(sort(hierr_object$penalty_ext))
    }
    if (loss %in% c("mse", "mae")) {
        opt_loss <- min(cv_mean, na.rm = TRUE)
        optIndex <- which(opt_loss == cv_mean, arr.ind = TRUE)
    } else {
        opt_loss <- max(cv_mean, na.rm = TRUE)
        optIndex <- which(opt_loss == cv_mean, arr.ind = TRUE)
    }

    if (is.null(dim(optIndex))) {
        opt_penalty <- hierr_object$penalty[optIndex[1]]
        opt_penalty_ext <- hierr_object$penalty_ext[optIndex[2]]
    } else {
        opt_penalty <- hierr_object$penalty[optIndex[1, 1]]
        opt_penalty_ext <- hierr_object$penalty_ext[optIndex[1, 2]]
    }

    cvfit <- list(cv_mean = cv_mean,
                  cv_sd = cv_sd,
                  opt_loss = opt_loss,
                  opt_penalty = opt_penalty,
                  opt_penalty_ext = opt_penalty_ext,
                  fit_train = hierr_object,
                  call = hierr_object$call)

    class(cvfit) <- c("cvhierr")
    return(cvfit)
}

