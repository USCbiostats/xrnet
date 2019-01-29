#' Predict function for "cvhierr" object
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
#' if p or pext are not in the original path(s) computed. See \code{\link{hierr}} for
#' more information on regularization object.
#' @param ... pass other arguments to hierr function (if needed)

#' @export
predict.cvhierr <- function(object,
                            newdata = NULL,
                            newdata_fixed = NULL,
                            p = NULL,
                            pext = NULL,
                            type = c("response", "coefficients"),
                            penalty = NULL,
                            ...)
{
    predict(object$fit_train,
            newdata = newdata,
            newdata_fixed = newdata_fixed,
            p = p,
            pext = pext,
            type = type,
            penalty = penalty,
            ...)
}
