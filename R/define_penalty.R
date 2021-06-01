#' Define regularization object for predictor and external data.
#'
#' @description Defines regularization for predictors and external data
#' variables in \code{\link{xrnet}} fitting. Use helper functions define_lasso,
#' define_ridge, or define_enet to specify a common penalty on x or external.
#'
#' @param penalty_type type of regularization. Default is 1 (Lasso).
#' Can supply either a scalar value or vector with length equal to the number of
#' variables the matrix.
#' \itemize{
#'    \item 0 = Ridge
#'    \item (0,1) = Elastic-Net
#'    \item 1 = Lasso / Quantile
#' }
#' @param quantile specifies quantile for quantile penalty. Default of 0.5
#' reduces to lasso (currently not implemented).
#' @param num_penalty number of penalty values to fit in grid. Default is 20.
#' @param penalty_ratio ratio between minimum and maximum penalty for x.
#' Default is 1e-04 if \eqn{n > p} and 0.01 if \eqn{n <= p}.
#' @param user_penalty user-defined vector of penalty values to use in penalty
#' path.
#' @param custom_multiplier variable-specific penalty multipliers to apply to
#' overall penalty. Default is 1 for all variables. 0 is no penalization.
#'
#' @return A list object with regularization settings that are used to define
#' the regularization for predictors or external data in \code{\link{xrnet}} and
#' \code{\link{tune_xrnet}}:
#' \item{penalty_type}{The penalty type, scalar with value in range [0, 1].}
#' \item{quantile}{Quantile for quantile penalty, 0.5 defaults to lasso
#' (not currently implemented).}
#' \item{num_penalty}{The number of penalty values in the penalty path.}
#' \item{penalty_ratio}{The ratio of the minimum penalty value compared to the
#' maximum penalty value.}
#' \item{user_penalty}{User-defined numeric vector of penalty values, NULL if
#' not provided by user.}
#' \item{custom_multiplier}{User-defined feature-specific penalty multipliers,
#' NULL if not provided by user.}
#'
#' @examples
#'
#' # define ridge penalty with penalty grid split into 30 values
#' my_penalty <- define_penalty(penalty_type = 0, num_penalty = 30)
#'
#' # define elastic net (0.5) penalty with user-defined penalty
#' my_custom_penalty <- define_penalty(
#'   penalty_type = 0.5, user_penalty = c(100, 50, 10, 1, 0.1)
#' )
#' @export
define_penalty <- function(penalty_type = 1,
                           quantile = 0.5,
                           num_penalty = 20,
                           penalty_ratio = NULL,
                           user_penalty = NULL,
                           custom_multiplier = NULL) {
  if (any(penalty_type < 0) || any(penalty_type > 1)) {
    stop("Invalid penalty type")
  } else {
    penalty_type <- as.double(penalty_type)
  }

  if (quantile < 0 || quantile > 1) {
    stop("Invalid value for quantile, must be between 0 and 1")
  } else {
    quantile <- as.double(quantile)
  }

  if (is.null(user_penalty)) {
    user_penalty <- as.double(0)
    num_penalty <- as.integer(num_penalty)
    if (!is.null(penalty_ratio)) {
      if (penalty_ratio <= 0 | penalty_ratio >= 1) {
        stop("penalty_ratio should be between 0 and 1")
      } else {
        penalty_ratio <- as.double(penalty_ratio)
      }
    }
  } else {
    penalty_ratio <- as.double(0)
    if (any(user_penalty < 0)) {
      stop("user_penalty can only contain non-negative values")
    }
    user_penalty <- as.double(rev(sort(user_penalty)))
    num_penalty <- as.integer(length(user_penalty))
  }

  if (!is.null(custom_multiplier) && any(custom_multiplier < 0)) {
    stop("custom_multiplier can only contain non-negative values")
  }

  penalty_obj <- list(
    penalty_type = penalty_type,
    quantile = quantile,
    num_penalty = num_penalty,
    penalty_ratio = penalty_ratio,
    user_penalty = user_penalty,
    custom_multiplier = custom_multiplier
  )
}

#' Define lasso regularization object for predictor and external data
#'
#' @description Helper function to define a lasso penalty regularization object.
#' See \code{define_penalty} for more details.
#'
#' @param num_penalty number of penalty values to fit in grid. Default is 20.
#' @param penalty_ratio ratio between minimum and maximum penalty for x.
#' Default is 1e-04 if \eqn{n > p} and 0.01 if \eqn{n <= p}.
#' @param user_penalty user-defined vector of penalty values to use in penalty
#' path.
#' @param custom_multiplier variable-specific penalty multipliers to apply to
#' overall penalty.
#' Default is 1 for all variables. 0 is no penalization.
#' @return A list object with regularization settings that are used to define
#' the regularization
#' for predictors or external data in \code{\link{xrnet}} and
#' \code{\link{tune_xrnet}}. The list
#' elements will match those returned by \code{\link{define_penalty}},
#' but with the penalty_type automatically set to 1.

#' @export
define_lasso <- function(num_penalty = 20,
                         penalty_ratio = NULL,
                         user_penalty = NULL,
                         custom_multiplier = NULL) {
  define_penalty(
    penalty_type = 1,
    quantile = 0.5,
    num_penalty = num_penalty,
    penalty_ratio = penalty_ratio,
    user_penalty = user_penalty,
    custom_multiplier = custom_multiplier
  )
}

#' Define ridge regularization object for predictor and external data
#'
#' @description Helper function to define a ridge penalty regularization object.
#' See \code{define_penalty} for more details.
#'
#' @param num_penalty number of penalty values to fit in grid. Default is 20.
#' @param penalty_ratio ratio between minimum and maximum penalty for x.
#' Default is 1e-04 if \eqn{n > p} and 0.01 if \eqn{n <= p}.
#' @param user_penalty user-defined vector of penalty values to use in penalty
#' path.
#' @param custom_multiplier variable-specific penalty multipliers to apply to
#' overall penalty. Default is 1 for all variables. 0 is no penalization.
#'
#' @return A list object with regularization settings that are used to define
#' the regularization for predictors or external data in \code{\link{xrnet}} and
#' \code{\link{tune_xrnet}}. The list elements will match those returned by
#' \code{\link{define_penalty}}, but with the penalty_type automatically set
#' to 0.

#' @export
define_ridge <- function(num_penalty = 20,
                         penalty_ratio = NULL,
                         user_penalty = NULL,
                         custom_multiplier = NULL) {
  define_penalty(
    penalty_type = 0,
    quantile = 0.5,
    num_penalty = num_penalty,
    penalty_ratio = penalty_ratio,
    user_penalty = user_penalty,
    custom_multiplier = custom_multiplier
  )
}

#' Define elastic net regularization object for predictor and external data
#'
#' @description Helper function to define a elastic net penalty regularization
#' object. See \code{define_penalty} for more details.
#'
#' @param en_param elastic net parameter, between 0 and 1
#' @param num_penalty number of penalty values to fit in grid. Default is 20.
#' @param penalty_ratio ratio between minimum and maximum penalty for x.
#' Default is 1e-04 if \eqn{n > p} and 0.01 if \eqn{n <= p}.
#' @param user_penalty user-defined vector of penalty values to use in penalty
#' path.
#' @param custom_multiplier variable-specific penalty multipliers to apply to
#' overall penalty. Default is 1 for all variables. 0 is no penalization.
#'
#' @return A list object with regularization settings that are used to define
#' the regularization for predictors or external data in \code{\link{xrnet}} and
#' \code{\link{tune_xrnet}}. The list elements will match those returned by
#' \code{\link{define_penalty}}, but with the penalty_type set to match the
#' value of \code{en_param}.

#' @export
define_enet <- function(en_param = 0.5,
                        num_penalty = 20,
                        penalty_ratio = NULL,
                        user_penalty = NULL,
                        custom_multiplier = NULL) {
  define_penalty(
    penalty_type = en_param,
    quantile = 0.5,
    num_penalty = num_penalty,
    penalty_ratio = penalty_ratio,
    user_penalty = user_penalty,
    custom_multiplier = custom_multiplier
  )
}
