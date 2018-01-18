#' k-fold cross-validation for hierarchical regularized regression
#'
#' @description k-fold cross-validation for hierarchical regularized regression \code{\link{hierr}}
#'
#' @param x predictor design matrix of dimension \eqn{n x p}
#' @param y outcome vector of length \eqn{n}
#' @param external external data design matrix of dimension \eqn{p x q}
#' @param family error distribution for outcome variable
#' @param penalty specifies regularization object for x and external. See \code{\link{definePenalty}} for more details.
#' @param weights optional vector of observation-specific weights. Default is 1 for all observations.
#' @param type.measure loss function for cross-validation. Options include:
#' \itemize{
#'    \item mse (Mean Squared Error)
#'    \item deviance
#'    \item mae (Mean Absolute Error)
#' }
#' @param nfolds number of folds for cross-validation. Default is 5.
#' @param parallel Use \code{foreach} function to fit folds in parallel if TRUE, must register cluster before using.
#' @param ... list of additional arguments to pass to function \code{\link{hierr}}.


#' @export
cvhierr <- function(x = x,
                    y = y,
                    external = ext,
                    family = c("gaussian", "binomial"),
                    penalty = definePenalty(),
                    weights = NULL,
                    type.measure = c("mse", "mae", "deviance", "class", "auc"),
                    nfolds = 5,
                    parallel = FALSE, ...)
{

    # Set measure used to assess model prediction performance
    if (missing(type.measure)) {
        type.measure <- "default"
    } else {
        type.measure <- match.arg(type.measure)
    }

    # Check family argument
    family <- match.arg(family)

    # Generate error function based on family/type.measure
    calc_error <- error_match(family = family, type.measure = type.measure)

    # Get arguments to cvhirr() function and filter for calls to fitting procedure
    hierr_call <- match.call(expand.dots = TRUE)
    cv_args <- match(c("type.measure", "nfolds", "foldid"), names(hierr_call), FALSE)

    if (any(cv_args)) {
        hierr_call <- hierr_call[-cv_args]
    }
    hierr_call[[1]] <- as.name("hierr")

    # Set sample size/ weights
    n <- nrow(x)
    if (is.null(weights)) {
        weights <- rep(1, n)
    }

    # Create hierr object
    hierr_object <- hierr(x = x, y = y, external = ext, family = family, weights = weights, penalty = penalty, ...)
    hierr_object$call <- hierr_call
    penalty_fixed <- definePenalty(penalty$penalty_type, penalty$penalty_type_ext, user_penalty = hierr_object$lam, user_penalty_ext = hierr_object$lam_ext)

    # Randomly sample observations into folds
    foldid <- sample(rep(seq(nfolds), length = n))

    # Vector to collect results of CV
    errormat <- matrix(NA, nrow = length(y), ncol = length(hierr_object$lam)*length(hierr_object$lam_ext))

    # Run k-fold CV
    if (parallel) {

    } else {
        for (k in seq_along(1:nfolds)) {

            # Split into test and train for k-th fold
            subset <- (foldid == k)
            if (is.vector(drop(y))) {
                y_train <- y[!subset]
            } else {
                y_train <- y[!subset, ]
            }
            x_train <- x[!subset, ]
            weights_train <- weights[!subset]

            # Fit model on k-th training fold
            fit_fold <- hierr(x = x_train, y = y_train, external = ext, weights = weights_train, penalty = penalty_fixed, ...)[c("beta0", "betas")]
            fit_fold$betas <- rbind(as.vector(t(fit_fold$beta0)), `dim<-`(aperm(fit_fold$betas, c(1, 3, 2)), c(dim(fit_fold$betas)[1], dim(fit_fold$betas)[2] * dim(fit_fold$betas)[3])))
            errormat[subset, ] <- calc_error(fit_fold$betas, y[subset], cbind(1, x[subset, ]), weights[subset])
        }
    }
    cv_mean <- apply(errormat, 2, stats::weighted.mean, w = weights)
    cv_sd <- sqrt(apply(sweep(errormat, 2L, cv_mean, FUN = "-")^2, 2, stats::weighted.mean, w = weights, na.rm = TRUE) / (length(y) - 1))
    cv_mean <- matrix(cv_mean, nrow = length(hierr_object$lam), byrow = TRUE)
    cv_sd <- matrix(cv_sd, nrow = length(hierr_object$lam), byrow = TRUE)
    rm(errormat)
    row.names(cv_mean) <- rev(sort(hierr_object$lam))
    colnames(cv_mean) <- rev(sort(hierr_object$lam_ext))
    row.names(cv_sd) <- rev(sort(hierr_object$lam))
    colnames(cv_sd) <- rev(sort(hierr_object$lam_ext))

    min_error <- min(cv_mean, na.rm = TRUE)
    optIndex <- which(min_error == cv_mean, arr.ind = TRUE)

    if (is.null(dim(optIndex))) {
        minl1 <- as.numeric(row.names(cv_mean)[optIndex[1]])
        minl2 <- as.numeric(colnames(cv_mean)[optIndex[2]])
    } else {
        minl1 <- as.numeric(row.names(cv_mean)[optIndex[1,1]])
        minl2 <- as.numeric(colnames(cv_mean)[optIndex[1,2]])
    }

    list(cv_mean = cv_mean,
         cv_sd = cv_sd,
         min_error = min_error,
         minl1 = minl1,
         minl2 = minl2,
         lambda1 = hierr_object$lam,
         lambda2 = hierr_object$lam_ext
        )
}

