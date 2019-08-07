#' Predict function for "xrnet" object
#'
#' @description Predict coefficients or response in new data using fitted model from an \code{\link{xrnet}} object
#'
#' @param object A \code{\link{xrnet}} object
#' @param newdata matrix with new values for penalized variables
#' @param newdata_fixed matrix with new values for unpenalized variables
#' @param p vector of penalty values to apply to predictor variables
#' @param pext vector of penalty values to apply to external data variables
#' @param type type of prediction to make using the xrnet model, options include
#' \itemize{
#'    \item coefficients
#'    \item response
#'    \item link (linear predictor)
#' }
#' @param penalty (optional) regularization object applied to original model object, only needed
#' if p or pext are not in the original path(s) computed. See \code{\link{define_penalty}} for
#' more information on regularization object.
#' @param ... pass other arguments to xrnet function (if needed)
#' @return The object returned based on the type object is as follows:
#' \itemize{
#'     \item coefficients: A list with the coefficient estimates for each penalty combination
#'     \item response: An array with the response predictions based on the data for each penalty combination
#'     \item link: An array with linear predictions based on the data for each penalty combination
#' }
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
#' \dontrun{
#' lambda1 <- fit_xrnet$penalty[10]
#' lambda2 <- fit_xrnet$penalty[10]
#'
#' coef_xrnet <- predict(
#'     fit_xrnet,
#'     p = lambda1,
#'     pext = lambda2,
#'     type = "coefficients"
#' )
#'
#' pred_xrnet <- predict(
#'     fit_xrnet,
#'     p = lambda1,
#'     pext = lambda2,
#'     newdata = x_linear,
#'     type = "response"
#' )
#' }

#' @export
predict.xrnet <- function(object,
                          newdata = NULL,
                          newdata_fixed = NULL,
                          p = NULL,
                          pext = NULL,
                          type = c("response", "coefficients", "link"),
                          penalty = NULL,
                          ...)
{

    if (missing(type)) {
        type <- "response"
    } else {
        type <- match.arg(type)
    }

    if (missing(newdata) && !match(type, c("coefficients"), FALSE)) {
        stop("Error: 'newdata' needs to be specified")
    }

    if (is.null(p)){
        stop("Error: p not specified")
    }

    if (!(all(p %in% object$penalty)) || !(all(pext %in% object$penalty_ext))) {

        if (is.null(penalty)) {
            stop("Error: Not all penalty values in path(s), please provide original regularization object")
        }

        if (!is.null(object$penalty_ext)) {
            if (is.null(pext)) {
                stop("Error: pext not specified")
            }
            penalty$user_penalty <- unique(rev(sort(c(object$penalty, p))))
            penalty$num_penalty <- length(penalty$user_penalty)
            penalty$user_penalty_ext <- unique(rev(sort(c(object$penalty_ext, pext))))
            penalty$num_penalty_ext <- length(penalty$user_penalty_ext)
        } else {
            penalty$user_penalty <- unique(rev(sort(c(object$penalty, p))))
            penalty$num_penalty <- length(penalty$user_penalty)

            if(!is.null(pext)) {
                warning(paste("Warning: No external data variables in model fit,
                              ignoring supplied external penalties pext = ", pext))
                pext <- NULL
            }
        }

        xrnet_call <- object$call
        xrnet_call[["penalty"]] <- as.name("penalty")
        add_args <- match.call(expand.dots = FALSE, call = xrnet_call)$...

        if (length(add_args)) {
            existing <- !is.na(match(names(add_args), names(xrnet_call)))
            for (arg in names(add_args[existing]))
                xrnet_call[[arg]] <- add_args[[arg]]
        }
        tryCatch(object <- eval(xrnet_call),
                 error = function(e) stop("Error: Unable to refit 'xrnet' object,
                                          please supply arguments used in original function call")
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
    betas <- object$betas[ , idxl1, idxl2, drop = F]
    gammas <- object$gammas[, idxl1, idxl2, drop = F]
    alpha0 <- object$alpha0[idxl1, idxl2, drop = F]
    alphas <- object$alphas[ , idxl1, idxl2, drop = F]

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
            if (!(typeof(newdata) %in% c("integer", "double"))) {
                stop("Error: newdata contains non-numeric values")
            }
            mattype_x <- 1
        }
        else if (is.big.matrix(newdata)) {
            if (!(bigmemory::describe(newdata)@description$type %in% c("integer", "double"))) {
                stop("Error: newdata contains non-numeric values")
            }
            mattype_x <- 2
        } else if ("dgCMatrix" %in% class(newdata)) {
            if (!(typeof(newdata@x) %in% c("integer", "double"))){
                stop("Error: newdata contains non-numeric values")
            }
            mattype_x <- 3
        } else {
            stop("Error: newdata must be a matrix, big.matrix, filebacked.big.matrix, or dgCMatrix")
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
            result <- aperm(
                a= array(
                      data = t(result),
                      c(length(pext), length(p),
                      dim(result)[1])
                ),
                perm = c(3, 2, 1)
            )
        } else {

        }
        return(drop(result))
    }
}

