#' Predict function for "cv_hierr" object
#'
#' @description Predict coefficients or response in new data
#'
#' @param object A \code{\link{cv_hierr}} object
#' @param newdata matrix with new values for penalized variables
#' @param newdata_fixed matrix with new values for unpenalized variables
#' @param p vector of penalty values to apply to predictor variables
#' @param pext vector of penalty values to apply to external data variables
#' @param type type of prediction to make using the hierr model
#' @param penalty regularization object applied to original model object, only needed
#' if p or pext are not in the original path(s) computed. See \code{\link{definePenalty}} for
#' more information on regularization object.
#' @param ... pass other arguments to hierr function (if needed)

#' @export
predict.cv_hierr <- function(object,
                            newdata = NULL,
                            newdata_fixed = NULL,
                            p = "opt",
                            pext = "opt",
                            type = c("response", "coefficients", "link"),
                            penalty = NULL,
                            ...)
{
    if (p == "opt")
        p <- object$opt_penalty
    if (pext == "opt")
        pext <- object$opt_penalty_ext

    predict(object$fit_train,
            newdata = newdata,
            newdata_fixed = newdata_fixed,
            p = p,
            pext = pext,
            type = type,
            penalty = penalty,
            ...)
}
