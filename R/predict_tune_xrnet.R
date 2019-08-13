#' Predict function for "tune_xrnet" object
#'
#' @description Predict coefficients or response in new data
#'
#' @param object A \code{\link{tune_xrnet}} object
#' @param newdata matrix with new values for penalized variables
#' @param newdata_fixed matrix with new values for unpenalized variables
#' @param p vector of penalty values to apply to predictor variables.
#' Default is optimal value in tune_xrnet object.
#' @param pext vector of penalty values to apply to external data variables.
#' Default is optimal value in tune_xrnet object.
#' @param type type of prediction to make using the xrnet model, options include:
#' \itemize{
#'    \item response
#'    \item link (linear predictor)
#'    \item coefficients
#' }
#' @param penalty (optional) regularization object applied to original model object, only needed
#' if p or pext are not in the original path(s) computed. See \code{\link{define_penalty}} for
#' more information on regularization object.
#' @param ... pass other arguments to xrnet function (if needed)

#' @export
predict.tune_xrnet <- function(object,
                               newdata = NULL,
                               newdata_fixed = NULL,
                               p = "opt",
                               pext = "opt",
                               type = c("response", "link", "coefficients"),
                               penalty = NULL,
                               ...)
{
    if (p == "opt")
        p <- object$opt_penalty
    if (pext == "opt")
        pext <- object$opt_penalty_ext

    predict(object$fitted_model,
            newdata = newdata,
            newdata_fixed = newdata_fixed,
            p = p,
            pext = pext,
            type = type,
            penalty = penalty,
            ...)
}
