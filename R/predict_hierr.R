#' Predict function for "hierr" object
#'
#' @description Predict coefficients or response in new data
#'
#' @param object A \code{\link{hierr}} object
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
#' @importFrom stats update
predict.hierr <- function(object,
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

    if (missing(newdata)) {
        if (!match(type, c("coefficients"), FALSE))
            stop("Error: 'newdata' needs to be specified")
    }

    if (!(all(p %in% object$penalty)) || !(all(pext %in% object$penalty_ext))) {

        if (is.null(penalty))
            stop("Error: Penalties not in original path, please provide original penalty object to refit model")

        if (!is.null(object$penalty_ext)) {
            if (is.null(p)) {
                stop("Error: p not specified")
            }
            if (is.null(pext)) {
                stop("Error: pext not specified")
            }
            penalty$user_penalty <- c(object$penalty, p)
            penalty$num_penalty <- length(penalty$user_penalty)
            penalty$user_penalty_ext <- c(object$penalty_ext, pext)
            penalty$num_penalty_ext <- length(penalty$user_penalty_ext)
        } else {
            if (is.null(p)) {
                stop("error: p not specified")
            }
            penalty$user_penalty <- c(object$penalty, p)
            penalty$num_penalty <- length(penalty$user_penalty)

            if(!is.null(pext)) {
                warning(paste("Warning: No external data variables in model fit, ignoring supplied external penalties pext = ", pext))
                pext <- NULL
            }
        }

        hierr_call <- object$call
        hierr_call[["penalty"]] <- as.name("penalty")

        tryCatch(object <- eval(hierr_call),
                 error = function(e) stop("Error: Unable to refit 'hierr' object, please supply arguments used in original function call")
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
            if (!(typeof(newdata) %in% c("integer", "double")))
                stop("Error: newdata contains non-numeric values")
            mattype_x <- 1
        }
        else if (is.big.matrix(newdata)) {
            if (!(bigmemory::describe(newdata)@description$type %in% c("integer", "double")))
                stop("Error: newdata contains non-numeric values")
            mattype_x <- 2
        } else if ("dgCMatrix" %in% class(newdata)) {
            if (!(typeof(newdata@x) %in% c("integer", "double")))
                stop("Error: newdata contains non-numeric values")
            mattype_x <- 3
        } else {
            stop("Error: newdata must be a matrix, big.matrix, filebacked.big.matrix, or dgCMatrix")
        }

        beta0 <- as.vector(beta0)
        betas <- `dim<-`(aperm(betas, c(1, 3, 2)), c(dim(betas)[1], dim(betas)[2] * dim(betas)[3]))
        if (!is.null(gammas)) {
            gammas <- `dim<-`(aperm(gammas, c(1, 3, 2)), c(dim(gammas)[1], dim(gammas)[2] * dim(gammas)[3]))
        } else {
            gammas <- matrix(vector("numeric", 0), 0, 0)
            newdata_fixed <- matrix(vector("numeric", 0), 0, 0)
        }

        if (mattype_x == 2)
            result <- computeResponseRcpp(newdata@address,
                                          mattype_x,
                                          newdata_fixed,
                                          beta0,
                                          betas,
                                          gammas,
                                          type,
                                          object$family)
        else
            result <- computeResponseRcpp(newdata,
                                          mattype_x,
                                          newdata_fixed,
                                          beta0,
                                          betas,
                                          gammas,
                                          type,
                                          object$family)

        if (length(pext) > 1) {
            result <- aperm(array(t(result), c(length(pext), length(p), dim(result)[1])), c(3, 2, 1))
        } else {

        }
        return(drop(result))
    }
}

