#' Predict function for "xrnet" object
#'
#' @description Extract coefficients or  predict response in new data using
#' fitted model from an \code{\link{xrnet}} object. Note that we currently only
#' support returning coefficient estimates that are in the original path(s).
#'
#' @param object A \code{\link{xrnet}} object
#' @param newdata matrix with new values for penalized variables
#' @param newdata_fixed matrix with new values for unpenalized variables
#' @param p vector of penalty values to apply to predictor variables
#' @param pext vector of penalty values to apply to external data variables
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
#' coef_xrnet <- predict(
#'   fit_xrnet,
#'   p = lambda1,
#'   pext = lambda2,
#'   type = "coefficients"
#' )
#'
#' pred_xrnet <- predict(
#'   fit_xrnet,
#'   p = lambda1,
#'   pext = lambda2,
#'   newdata = x_linear,
#'   type = "response"
#' )
#' @export
predict.xrnet <- function(object,
                          newdata = NULL,
                          newdata_fixed = NULL,
                          p = NULL,
                          pext = NULL,
                          type = c("response", "link", "coefficients"),
                          ...) {
  if (missing(type)) {
    type <- "response"
  } else {
    type <- match.arg(type)
  }

  if (missing(newdata) && !match(type, c("coefficients"), FALSE)) {
    stop("newdata needs to be specified")
  }

  if (is.null(p)) {
    stop("p not specified")
  }

  if (!is.null(object$penalty_ext) && is.null(pext)) {
    stop("pext not specified")
  }

  if (!(all(p %in% object$penalty)) || !(all(pext %in% object$penalty_ext))) {
    stop(
      "Not all penalty values in path(s),
      please refit xrnet() model with desired penalty values"
    )
  }

  p <- rev(sort(p))
  idxl1 <- which(object$penalty %in% p)
  if (!is.null(object$penalty_ext)) {
    pext <- rev(sort(pext))
    idxl2 <- which(object$penalty_ext %in% pext)
  } else {
    idxl2 <- 1
  }

  beta0 <- object$beta0[idxl1, idxl2, drop = F]
  betas <- object$betas[, idxl1, idxl2, drop = F]
  gammas <- object$gammas[, idxl1, idxl2, drop = F]
  alpha0 <- object$alpha0[idxl1, idxl2, drop = F]
  alphas <- object$alphas[, idxl1, idxl2, drop = F]

  if (type == "coefficients") {
    return(list(
      beta0 = beta0,
      betas = betas,
      gammas = gammas,
      alpha0 = alpha0,
      alphas = alphas,
      penalty = p,
      penalty_ext = pext
    ))
  }

  if (type %in% c("link", "response")) {
    if (is(newdata, "matrix")) {
      if (typeof(newdata) != "double") {
        stop("newdata must be of type double")
      }
      mattype_x <- 1
    }
    else if (is.big.matrix(newdata)) {
      if (bigmemory::describe(newdata)@description$type != "double") {
        stop("newdata must be of type double")
      }
      mattype_x <- 2
    } else if ("dgCMatrix" %in% class(newdata)) {
      if (typeof(newdata@x) != "double") {
        stop("newdata must be of type double")
      }
      mattype_x <- 3
    } else {
      stop(
        "newdata must be a matrix, big.matrix,
        filebacked.big.matrix, or dgCMatrix"
      )
    }

    beta0 <- as.vector(beta0)
    betas <- `dim<-`(
      aperm(betas, c(1, 3, 2)),
      c(dim(betas)[1], dim(betas)[2] * dim(betas)[3])
    )
    if (!is.null(gammas)) {
      gammas <- `dim<-`(
        aperm(gammas, c(1, 3, 2)),
        c(dim(gammas)[1], dim(gammas)[2] * dim(gammas)[3])
      )
    } else {
      gammas <- matrix(vector("numeric", 0), 0, 0)
      newdata_fixed <- matrix(vector("numeric", 0), 0, 0)
    }

    result <- computeResponseRcpp(
      newdata,
      mattype_x,
      newdata_fixed,
      beta0,
      betas,
      gammas,
      type,
      object$family
    )

    if (length(pext) > 1) {
      dim(result) <- c(NROW(result), length(pext), length(p))
      result <- aperm(result, c(1, 3, 2))
    }
    return(drop(result))
  }
}
