#' k-fold cross-validation for hierarchical regularized regression
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#'
#' @description k-fold cross-validation for hierarchical regularized regression \code{\link{hierr}}
#'
#' @param x predictor design matrix of dimension \eqn{n x p}
#' @param y outcome vector of length \eqn{n}
#' @param external (optional) external data design matrix of dimension \eqn{p x q}
#' @param unpen (optional) unpenalized predictor design matrix
#' @param family error distribution for outcome variable
#' @param penalty specifies regularization object for x and external. See \code{\link{definePenalty}} for more details.
#' @param weights optional vector of observation-specific weights. Default is 1 for all observations.
#' @param loss loss function for cross-validation. Options include:
#' \itemize{
#'    \item mse (Mean Squared Error)
#'    \item mae (Mean Absolute Error)
#'    \item auc (Area under the curve)
#' }
#' @param nfolds number of folds for cross-validation. Default is 5.
#' @param foldid (optional) vector that identifies user-specified fold for each observation. If NULL, folds are automatically generated.
#' @param parallel use \code{foreach} function to fit folds in parallel if TRUE, must register cluster (\code{doParallel}) before using.
#' @param ... list of additional arguments to pass to function \code{\link{hierr}}.

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

    # Set measure used to assess model prediction performance
    if (missing(loss)) {
        loss <- "default"
    } else {
        loss <- match.arg(loss)
    }

    # Check family argument
    family <- match.arg(family)

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
    }
    if (is.null(external)) {
        external <- matrix(vector("numeric", 0), 0, 0)
        nc_ext <- as.integer(0)
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
                                          quantiles = c(penalty_fold$tau, penalty_fold$tau_ext),
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
                                          quantiles = c(penalty_fold$tau, penalty_fold$tau_ext),
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
                                          quantiles = c(penalty_fold$tau, penalty_fold$tau_ext),
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
        optl1 <- hierr_object$penalty[optIndex[1]]
        optl2 <- hierr_object$penalty_ext[optIndex[2]]
    } else {
        optl1 <- hierr_object$penalty[optIndex[1, 1]]
        optl2 <- hierr_object$penalty_ext[optIndex[1, 2]]
    }

    cvfit <- list(cv_mean = cv_mean,
                  cv_sd = cv_sd,
                  opt_loss = opt_loss,
                  optl1 = optl1,
                  optl2 = optl2,
                  penalty = hierr_object$penalty,
                  penalty_ext = hierr_object$penalty_ext,
                  hierr_fit = hierr_object,
                  call = hierr_object$call)

    class(cvfit) <- c("cvhierr")
    return(cvfit)
}

