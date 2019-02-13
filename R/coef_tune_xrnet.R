#' Get coefficient estimates from "tune_xrnet" model object
#'
#' @description Returns coefficents from 'xrnet' model or refits model to estimate
#' coefficients not computed in path(s).
#'
#' @param object A \code{\link{tune_xrnet}} object
#' @param p vector of penalty values to apply to predictor variables.
#' Default is optimal value in tune_xrnet object.
#' @param pext vector of penalty values to apply to external data variables.
#' Default is optimal value in tune_xrnet object.
#' @param penalty (optional) regularization object applied to original model object, only needed
#' if p or pext are not in the original path(s) computed. See \code{\link{define_penalty}} for
#' more information on regularization object.
#' @param ... pass other arguments to xrnet function (if needed)

#' @export
coef.tune_xrnet <- function(object,
                          p = "opt",
                          pext = "opt",
                          penalty = NULL,
                          ...) {

    if (p == "opt")
        p <- object$opt_penalty
    if (pext == "opt")
        pext <- object$opt_penalty_ext

    predict(object$fitted_model,
            newdata = NULL,
            newdata_fixed = NULL,
            p = p,
            pext = pext,
            type = "coefficients",
            penalty = penalty,
            ...)

}
