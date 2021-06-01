#' Get coefficient estimates from "tune_xrnet" model object.
#'
#' @description Returns coefficients from 'xrnet' model. Note that we currently
#' only support returning coefficient estimates that are in the original
#' path(s).
#'
#' @param object A \code{\link{tune_xrnet}} object.
#' @param p vector of penalty values to apply to predictor variables.
#' Default is optimal value in tune_xrnet object.
#' @param pext vector of penalty values to apply to external data variables.
#' Default is optimal value in tune_xrnet object.
#' @param ... pass other arguments to xrnet function (if needed).
#'
#' @return A list with coefficient estimates at each of the requested penalty
#' combinations.
#' \item{beta0}{matrix of first-level intercepts indexed by penalty values, NULL
#' if no first-level intercept in original model fit.}
#' \item{betas}{3-dimensional array of first-level penalized coefficients
#' indexed by penalty values.}
#' \item{gammas}{3-dimensional array of first-level non-penalized coefficients
#' indexed by penalty values, NULL if unpen NULL in original model fit.}
#' \item{alpha0}{matrix of second-level intercepts indexed by penalty values,
#' NULL if no second-level intercept in original model fit.}
#' \item{alphas}{3-dimensional array of second-level external data coefficients
#' indexed by penalty values, NULL if external NULL in original model fit.}
#'
#' @examples
#' ## Cross validation of hierarchical linear regression model
#' data(GaussianExample)
#'
#' ## 5-fold cross validation
#' cv_xrnet <- tune_xrnet(
#'   x = x_linear,
#'   y = y_linear,
#'   external = ext_linear,
#'   family = "gaussian",
#'   control = xrnet_control(tolerance = 1e-6)
#' )
#'
#' ## Get coefficient estimates at optimal penalty combination
#' coef_opt <- coef(cv_xrnet)
#' @export
coef.tune_xrnet <- function(object,
                            p = "opt",
                            pext = "opt",
                            ...) {
  if (p == "opt") {
    p <- object$opt_penalty
  }
  if (pext == "opt") {
    pext <- object$opt_penalty_ext
  }

  predict(
    object$fitted_model,
    newdata = NULL,
    newdata_fixed = NULL,
    p = p,
    pext = pext,
    type = "coefficients",
    ...
  )
}
