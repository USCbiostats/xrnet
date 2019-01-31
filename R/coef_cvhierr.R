#' Get coefficient estimates from "cv_hierr" model object
#'
#' @description Returns coefficents from 'hierr' model or refits model to estimate
#' coefficients not computed in path(s).
#'
#' @param object A \code{\link{cv_hierr}} object
#' @param p vector of penalty values to apply to predictor variables.
#' Default is optimal value in cv_hierr object.
#' @param pext vector of penalty values to apply to external data variables.
#' Default is optimal value in cv_hierr object.
#' @param penalty (optional) regularization object applied to original model object, only needed
#' if p or pext are not in the original path(s) computed. See \code{\link{define_penalty}} for
#' more information on regularization object.
#' @param ... pass other arguments to hierr function (if needed)

#' @export
coef.cv_hierr <- function(object,
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
