#' Predict function for "tune_xrnet" object
#'
#' @description Extract coefficients or predict response in new data using
#' fitted model from a \code{\link{tune_xrnet}} object. Note that we currently
#' only support returning results that are in the original path(s).
#'
#' @param object A \code{\link{tune_xrnet}} object
#' @param newdata matrix with new values for penalized variables
#' @param newdata_fixed matrix with new values for unpenalized variables
#' @param p vector of penalty values to apply to predictor variables.
#' Default is optimal value in tune_xrnet object.
#' @param pext vector of penalty values to apply to external data variables.
#' Default is optimal value in tune_xrnet object.
#' @param type type of prediction to make using the xrnet model, options
#' include:
#' \itemize{
#'    \item response
#'    \item link (linear predictor)
#'    \item coefficients
#' }
#' @param ... pass other arguments to xrnet function (if needed)
#'
#' @return The object returned is based on the value of type as follows:
#' \itemize{
#'     \item response: An array with the response predictions based on the data
#'     for each penalty combination
#'     \item link: An array with linear predictions based on the data for each
#'     penalty combination
#'     \item coefficients: A list with the coefficient estimates for each
#'     penalty combination. See \code{\link{coef.xrnet}}.
#' }
#'
#' @examples
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
#' ## Get coefficients and predictions at optimal penalty combination
#' coef_xrnet <- predict(cv_xrnet, type = "coefficients")
#' pred_xrnet <- predict(cv_xrnet, newdata = x_linear, type = "response")
#' @export
predict.tune_xrnet <- function(object,
                               newdata = NULL,
                               newdata_fixed = NULL,
                               p = "opt",
                               pext = "opt",
                               type = c("response", "link", "coefficients"),
                               ...) {
  if (p == "opt") {
    p <- object$opt_penalty
  }
  if (pext == "opt") {
    pext <- object$opt_penalty_ext
  }

  predict(object$fitted_model,
    newdata = newdata,
    newdata_fixed = newdata_fixed,
    p = p,
    pext = pext,
    type = type,
    ...
  )
}
