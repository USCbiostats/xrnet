#' Get coefficient estimates from "xrnet" model object
#'
#' @description Returns coefficents from 'xrnet' model or refit model to estimate
#' coefficients not computed in original path(s) of model object.
#'
#' @param object A \code{\link{xrnet}} object
#' @param p vector of penalty values to apply to predictor variables.
#' @param pext vector of penalty values to apply to external data variables.
#' @param penalty (optional) regularization object applied to original model object, only needed
#' if p or pext are not in the original path(s) computed. See \code{\link{define_penalty}} for
#' more information on regularization object.
#' @param ... pass other arguments to xrnet function (if needed)
#' @return A list with coefficient estimates at each of the requested penalty combinations
#' @examples
#' data(GaussianExample)
#'
#' fit_xrnet <- xrnet(
#'     x = x_linear,
#'     y = y_linear,
#'     external = ext_linear,
#'     family = "gaussian"
#' )
#'
#' lambda1 <- fit_xrnet$penalty[10]
#' lambda2 <- fit_xrnet$penalty[10]
#' \dontrun{
#' coef_xrnet <- coef(
#'     fit_xrnet,
#'     p = lambda1,
#'     pext = lambda2,
#' )
#' }
#'

#' @export
coef.xrnet <- function(object,
                       p = NULL,
                       pext = NULL,
                       penalty = NULL,
                       ...) {
    predict(
        object,
        newdata = NULL,
        newdata_fixed = NULL,
        p = p,
        pext = pext,
        type = "coefficients",
        penalty = penalty,
        ...
    )
}
