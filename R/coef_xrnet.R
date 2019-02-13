#' Get coefficient estimates from "xrnet" model object
#'
#' @description Returns coefficents from 'xrnet' model or refits model to estimate
#' coefficients not computed in path(s).
#'
#' @param object A \code{\link{xrnet}} object
#' @param p vector of penalty values to apply to predictor variables.
#' @param pext vector of penalty values to apply to external data variables.
#' @param penalty (optional) regularization object applied to original model object, only needed
#' if p or pext are not in the original path(s) computed. See \code{\link{define_penalty}} for
#' more information on regularization object.
#' @param ... pass other arguments to xrnet function (if needed)

#' @export
coef.xrnet <- function(object,
                       p = NULL,
                       pext = NULL,
                       penalty = NULL,
                       ...) {
    predict(object,
            newdata = NULL,
            newdata_fixed = NULL,
            p = p,
            pext = pext,
            type = "coefficients",
            penalty = penalty,
            ...)

}
