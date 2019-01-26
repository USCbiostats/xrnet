#' Predict function for "cvhierr" object
#'
#' @description Predict coefficients or response in new data
#'
#' @param object A \code{\link{cv_hierr}} object
#' @param newdata matrix with new X values
#' @param p vector of penalty values to apply to predictor variables
#' @param pext vector of penalty values to apply to external data variables
#' @param type type of prediction to make using the hierr model
#' @param ... pass other arguments to hierr function (if needed)

#' @export
predict.cvhierr <- function(object,
                            newdata = NULL,
                            p = NULL,
                            pext = NULL,
                            type = c("response", "coefficients"),
                            ...)
{
    predict(object$hierr_fit, newdata = newdata, p = p, pext = pext, type = type, ...)
}
