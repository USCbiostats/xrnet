#' Define regularization object for predictor and external data
#'
#' @description Defines regularization terms for predictor and external data in \code{\link{xrnet}} fitting.
#'
#' @param penalty_type type of regularization for x. Default is 0 (Ridge). Can supply either a scalar value or vector with length equal to the number of variables in x.
#' \itemize{
#'    \item 0 = Ridge
#'    \item (0,1) = Elastic-Net
#'    \item 1 = Lasso / Quantile
#' }
#' @param quantile specifies quantile for predictors. Default of 0.5 reduces to lasso.
#' @param penalty_type_ext type of regularization for external data. See penalty_type for options. Default is 1 (lasso). Can supply either a scalar value or vector with length equal to the number of variables in external.
#' @param quantile_ext specifies quantile for external data. Default of 0.5 reduces to lasso.
#' @param num_penalty number of penalty values to fit in grid for x. Default is 20.
#' @param num_penalty_ext number of penalty values to fit in grid for external data. Default is 20.
#' @param penalty_ratio ratio between minimum and maximum penalty for x. Default is 1e-04 if \eqn{n > p} and 0.01 if \eqn{n <= p}.
#' @param penalty_ratio_ext ratio between minimum and maximum penalty for external data. Default is 1e-04 if \eqn{p > q} and 0.01 if \eqn{p <= q}.
#' @param user_penalty user-defined vector of penalty values to fit for x.
#' @param user_penalty_ext user-defined vector of penalty values to fit for external data.
#' @param custom_multiplier variable-specific penalty multipliers for x. Default is 1 for all variables. 0 is no penalization.
#' @param custom_multiplier_ext variable-specific penalty multipliers for external data. Default is 1 for all variables. 0 is no penalization.

#' @export
define_penalty <- function(penalty_type = 0,
                           quantile = 0.5,
                           penalty_type_ext = 1,
                           quantile_ext = 0.5,
                           num_penalty = 20,
                           num_penalty_ext = 20,
                           penalty_ratio = NULL,
                           penalty_ratio_ext = NULL,
                           user_penalty = NULL,
                           user_penalty_ext = NULL,
                           custom_multiplier = NULL,
                           custom_multiplier_ext = NULL) {

    if (any(penalty_type < 0) || any(penalty_type > 1)) {
        stop("Error: Invalid penalty type for x")
    } else {
        penalty_type <- as.double(penalty_type)
    }

    if (quantile < 0 || quantile > 1) {
        stop("Error: invalid value for quantile, must be between 0 and 1")
    } else {
        quantile <- as.double(quantile)
    }

    if (any(penalty_type_ext < 0) || any(penalty_type_ext > 1)) {
        stop("Error: Invalid penalty type for external")
    } else {
        penalty_type_ext <- as.double(penalty_type_ext)
    }

    if (quantile_ext < 0 || quantile_ext > 1) {
        stop("Error: invalid value for quantile_ext, must be between 0 and 1")
    } else {
        quantile_ext <- as.double(quantile_ext)
    }

    if (is.null(user_penalty)) {
        user_penalty <- as.double(0)
        num_penalty <- as.integer(num_penalty)
        if (!is.null(penalty_ratio)) {
            if (penalty_ratio <= 0 | penalty_ratio >= 1) {
                stop("Error: penalty_ratio should be between 0 and 1")
            } else {
                penalty_ratio <- as.double(penalty_ratio)
            }
        }
    } else {
        penalty_ratio <- as.double(0)
        if (any(user_penalty < 0)) {
            stop("Error: user_penalty can only contain non-negative values")
        }
        user_penalty <- as.double(rev(sort(user_penalty)))
        num_penalty <- as.integer(length(user_penalty))
    }

    if (is.null(user_penalty_ext)) {
        user_penalty_ext <- as.double(0)
        num_penalty_ext <- as.integer(num_penalty_ext)
        if (!is.null(penalty_ratio_ext)) {
            if (penalty_ratio_ext <=0 | penalty_ratio_ext >= 1) {
                stop("Error: penalty_ratio_ext should be between 0 and 1")
            } else {
                penalty_ratio_ext <- as.double(penalty_ratio_ext)
            }
        }
    } else {
        penalty_ratio_ext <- as.double(0)
        if (any(user_penalty_ext < 0)) {
            stop("Error: user_penalty_ext can only contain non-negative values")
        }
        user_penalty_ext <- as.double(rev(sort(user_penalty_ext)))
        num_penalty_ext <- as.integer(length(user_penalty_ext))
    }

    if (!is.null(custom_multiplier) && any(custom_multiplier < 0)) {
        stop("Error: custom_multiplier can only contain non-negative values")
    }

    if (!is.null(custom_multiplier_ext) && any(custom_multiplier_ext < 0)) {
        stop("Error: custom_multiplier_ext can only contain non-negative values")
    }

    penalty_obj <- list(penalty_type = penalty_type,
                        quantile = quantile,
                        penalty_type_ext = penalty_type_ext,
                        quantile_ext = quantile_ext,
                        num_penalty = num_penalty,
                        num_penalty_ext = num_penalty_ext,
                        penalty_ratio = penalty_ratio,
                        penalty_ratio_ext = penalty_ratio_ext,
                        user_penalty = user_penalty,
                        user_penalty_ext = user_penalty_ext,
                        custom_multiplier = custom_multiplier,
                        custom_multiplier_ext = custom_multiplier_ext)
}
