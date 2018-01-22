#' Predict function for "hierr" object
#'
#' @description Predict coefficients or response in new data
#'
#' @param object A \code{\link{hierr}} object
#' @param newdata matrix with new X values
#' @param p vector of penalty values to apply to predictor variables
#' @param pext vector of penalty values to apply to external data variables
#' @param type type of prediction to make using the hierr model
#' @param ... pass other arguments to hierr function (if needed)

#' @export
#' @importFrom stats update
predict.hierr <- function(object,
                         newdata = NULL,
                         p = NULL,
                         pext = NULL,
                         type = c("response", "coefficients"),
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

    if (!(all(p %in% object$penalty > 0)) | !(all(pext %in% object$penalty_ext > 0))) {
        tryCatch(object <- update(object, penalty = definePenalty(penalty_type = object$penalty_type,
                                                                  penalty_type_ext = object$penalty_type_ext,
                                                                  user_penalty = c(object$penalty, p),
                                                                  user_penalty_ext = c(object$penalty, pext)), ...),
                 error = function(e) stop("Error: Unable to refit 'hierr' object, please supply arguments used in original function call")
        )
    }

    idxl1 <- which(object$penalty %in% p)
    idxl2 <- which(object$penalty_ext %in% pext)

    beta0 <- object$beta0[idxl1, idxl2, drop = F]
    betas <- object$betas[ , idxl1, idxl2, drop = F]
    alpha0 <- object$alpha0[idxl1, idxl2, drop = F]
    alphas <- object$alphas[ , idxl1, idxl2, drop = F]
    p <- rev(sort(p))
    pext <- rev(sort(pext))

    if (type == "coefficients") {
        return(list(
            beta0 = beta0,
            betas = betas,
            alpha0 = alpha0,
            alphas = alphas,
            penalty = p,
            penalty_ext = pext,
            penalty_type = object$penalty_type,
            penalty_type_ext = object$penalty_type_ext
        ))
    }

    if (type == "response") {
        betas <- rbind(as.vector(t(beta0)), `dim<-`(aperm(betas, c(1, 3, 2)), c(dim(betas)[1], dim(betas)[2] * dim(betas)[3])))
        result <- cbind(1, newdata) %*% betas
        result <- aperm(array(t(result), c(length(pext), length(p), dim(result)[1])), c(3, 2, 1))
        return(drop(result))
    }
}

