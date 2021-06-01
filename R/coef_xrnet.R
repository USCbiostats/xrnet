#' Get coefficient estimates from "xrnet" model object.
#'
#' @description Returns coefficients from 'xrnet' model. Note that we currently
#' only support returning coefficient estimates that are in the original
#' path(s).
#'
#' @param object A \code{\link{xrnet}} object.
#' @param p vector of penalty values to apply to predictor variables.
#' @param pext vector of penalty values to apply to external data variables.
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
#' data(GaussianExample)
#'
#' fit_xrnet <- xrnet(
#'   x = x_linear,
#'   y = y_linear,
#'   external = ext_linear,
#'   family = "gaussian"
#' )
#'
#' lambda1 <- fit_xrnet$penalty[10]
#' lambda2 <- fit_xrnet$penalty_ext[10]
#'
#' coef_xrnet <- coef(
#'   fit_xrnet,
#'   p = lambda1,
#'   pext = lambda2,
#' )
#' @export
coef.xrnet <- function(object,
                       p = NULL,
                       pext = NULL,
                       ...) {
  predict(
    object,
    newdata = NULL,
    newdata_fixed = NULL,
    p = p,
    pext = pext,
    type = "coefficients",
    ...
  )
}
