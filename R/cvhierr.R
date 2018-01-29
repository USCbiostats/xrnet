#' k-fold cross-validation for hierarchical regularized regression
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
#' @param type.measure loss function for cross-validation. Options include:
#' \itemize{
#'    \item mse (Mean Squared Error)
#'    \item deviance
#'    \item mae (Mean Absolute Error)
#' }
#' @param nfolds number of folds for cross-validation. Default is 5.
#' @param parallel Use \code{foreach} function to fit folds in parallel if TRUE, must register cluster (\code{doParallel}) before using.
#' @param ... list of additional arguments to pass to function \code{\link{hierr}}.


#' @export
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
cvhierr <- function(x,
                    y,
                    external = NULL,
                    unpen = NULL,
                    family = c("gaussian"),
                    penalty = definePenalty(0, 1),
                    weights = NULL,
                    type.measure = c("mse", "mae", "deviance"),
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

    # Set sample size / weights
    n <- length(y)
    if (is.null(weights)) {
        weights <- rep(1, n)
    }

    # Create hierr object
    hierr_object <- hierr(x = x, y = y, external = external, unpen = unpen, family = family, weights = weights, penalty = penalty, ...)
    hierr_object$call <- hierr_call
    penalty_fixed <- definePenalty(penalty$penalty_type,
                                   penalty$penalty_type_ext,
                                   user_penalty = hierr_object$penalty,
                                   user_penalty_ext = hierr_object$penalty_ext)

    num_pen <- length(hierr_object$penalty)
    if (!is.null(external)) {
        num_pen_ext <- length(hierr_object$penalty_ext)
    } else {
        num_pen_ext <- 1
    }

    # Randomly sample observations into folds
    foldid <- sample(rep(seq(nfolds), length = n))

    # Vector to collect results of CV
    errormat <- matrix(NA, nrow = n, ncol = num_pen * num_pen_ext)

    # Run k-fold CV
    if (parallel) {
        cvout <- foreach(k = seq(nfolds), .packages = c("hierr")) %dopar% {
            subset <- (foldid == k)
            if (is.vector(drop(y))) {
                y_train <- y[!subset]
            } else {
                y_train <- y[!subset, ]
            }
            x_train <- x[!subset, ]
            if (!is.null(unpen)) {
                unpen_train <- unpen[!subset, ]
            } else {
                unpen_train <- NULL
            }
            weights_train <- weights[!subset]

            # Fit model on k-th training fold
            hierr(x = x_train, y = y_train, external = external, unpen = unpen_train, weights = weights_train, family = family, penalty = penalty_fixed, ...)[c("beta0", "betas")]
        }
        for (k in 1:nfolds) {
            subset <- (foldid == k)
            betas <- rbind(as.vector(t(cvout[[k]]$beta0)), `dim<-`(aperm(cvout[[k]]$betas, c(1, 3, 2)), c(dim(cvout[[k]]$betas)[1], dim(cvout[[k]]$betas)[2] * dim(cvout[[k]]$betas)[3])))
            errormat[subset, ] <- calc_error(betas, y[subset], cbind(1, x[subset, ]), weights[subset])

        }

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
            if (!is.null(unpen)) {
                unpen_train <- unpen[!subset, ]
            } else {
                unpen_train <- NULL
            }
            weights_train <- weights[!subset]

            # Fit model on k-th training fold
            fit_fold <- hierr(x = x_train, y = y_train, external = external, unpen = unpen_train, weights = weights_train, family = family, penalty = penalty_fixed, ...)[c("beta0", "betas", "gammas")]
            if (!is.null(fit_fold$gammas)) {
                betas <- rbind(as.vector(t(fit_fold$beta0)),
                               `dim<-`(aperm(fit_fold$betas, c(1, 3, 2)), c(dim(fit_fold$betas)[1], dim(fit_fold$betas)[2] * dim(fit_fold$betas)[3])),
                               `dim<-`(aperm(fit_fold$gammas, c(1, 3, 2)), c(dim(fit_fold$gammas)[1], dim(fit_fold$gammas)[2] * dim(fit_fold$gammas)[3])))
            } else {
                betas <- rbind(as.vector(t(fit_fold$beta0)),
                               `dim<-`(aperm(fit_fold$betas, c(1, 3, 2)), c(dim(fit_fold$betas)[1], dim(fit_fold$betas)[2] * dim(fit_fold$betas)[3])))
            }
            errormat[subset, ] <- calc_error(betas, y[subset], cbind(1, x[subset, ], unpen[subset, ]), weights[subset])
        }
    }
    cv_mean <- apply(errormat, 2, stats::weighted.mean, w = weights)
    cv_sd <- sqrt(apply(sweep(errormat, 2L, cv_mean, FUN = "-")^2, 2, stats::weighted.mean, w = weights, na.rm = TRUE) / (n - 1))
    cv_mean <- matrix(cv_mean, nrow = num_pen, byrow = TRUE)
    cv_sd <- matrix(cv_sd, nrow = num_pen, byrow = TRUE)
    row.names(cv_mean) <- rev(sort(hierr_object$penalty))
    colnames(cv_mean) <- rev(sort(hierr_object$penalty_ext))
    row.names(cv_sd) <- rev(sort(hierr_object$penalty))
    colnames(cv_sd) <- rev(sort(hierr_object$penalty_ext))

    min_error <- min(cv_mean, na.rm = TRUE)
    optIndex <- which(min_error == cv_mean, arr.ind = TRUE)

    if (is.null(dim(optIndex))) {
        minl1 <- as.numeric(row.names(cv_mean)[optIndex[1]])
        minl2 <- as.numeric(colnames(cv_mean)[optIndex[2]])
    } else {
        minl1 <- as.numeric(row.names(cv_mean)[optIndex[1,1]])
        minl2 <- as.numeric(colnames(cv_mean)[optIndex[1,2]])
    }

    cvfit <- list(cv_mean = cv_mean,
                  cv_sd = cv_sd,
                  min_error = min_error,
                  minl1 = minl1,
                  minl2 = minl2,
                  penalty = hierr_object$penalty,
                  penalty_ext = hierr_object$penalty_ext,
                  call = hierr_object$call)

    class(cvfit) <- c("cvhierr", "hierr")
    return(cvfit)
}

