#' Plot k-fold cross-validation error grid
#'
#' @description Generates plots to visualize the mean cross-validation error.
#' If no external data was used in the model fit, a plot of the cross-validated
#' error with standard error bars is generated for all penalty values. If
#' external data was used in the model fit, a contour plot of the
#' cross-validated errors is created. Error curves can also be generated for a
#' fixed value of the primary penalty on x (p) or the external penalty (pext)
#' when external data is used.
#'
#' @param x A tune_xrnet class object
#' @param p (optional) penalty value for x (for generating an error curve across
#' external penalties). Use value "opt" to use the optimal penalty value.
#' @param pext (optional) penalty value for external (for generating an error
#' curve across primary penalties). Use value "opt" to use the optimal penalty
#' value.
#' @param ... Additional graphics parameters
#'
#' @return None
#'
#' @details The parameter values p and pext can be used to generate profiled
#' error curves by fixing either the penalty on x or the penalty on external to
#' a fixed value. You cannot specify both at the same time as this would only
#' return a single point.
#'
#' @examples
#'
#' ## load example data
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
#' ## contour plot of cross-validated error
#' plot(cv_xrnet)
#'
#' ## error curve of external penalties at optimal penalty value
#' plot(cv_xrnet, p = "opt")
#' @export
#' @importFrom graphics filled.contour axis points
#' @importFrom grDevices colorRampPalette

plot.tune_xrnet <- function(x, p = NULL, pext = NULL, ...) {
  if (is.null(x$fitted_model$alphas) || !is.null(p) || !is.null(pext)) {
    if (is.null(x$fitted_model$alphas)) {
      xval <- log(as.numeric(rownames(x$cv_mean)))
      cverr <- x$cv_mean[, 1]
      cvsd <- x$cv_sd[, 1]
      xlab <- "log(Penalty)"
      xopt_val <- log(x$opt_penalty)
    } else {
      if (!is.null(p) && !is.null(pext)) {
        stop(
          "Please only specify either penalty or penalty_ext,
          cannot specify both at the same time"
        )
      } else if (!is.null(p)) {
        if (p == "opt") {
          p <- x$opt_penalty
        }
        p_idx <- match(p, x$fitted_model$penalty)
        if (is.na(p_idx)) {
          stop("The penalty value 'p' is not in the fitted model")
        }
        xval <- log(as.numeric(colnames(x$cv_mean)))
        cverr <- x$cv_mean[p_idx, ]
        cvsd <- x$cv_sd[p_idx, ]
        xlab <- "log(External Penalty)"
        xopt_val <- log(x$opt_penalty_ext)
      } else {
        if (pext == "opt") {
          pext <- x$opt_penalty_ext
        }
        pext_idx <- match(pext, x$fitted_model$penalty_ext)
        if (is.na(pext_idx)) {
          stop("The penalty value 'p' is not in the fitted model")
        }
        xval <- log(as.numeric(rownames(x$cv_mean)))
        cverr <- x$cv_mean[, pext_idx]
        cvsd <- x$cv_sd[, pext_idx]
        xlab <- "log(Penalty)"
        xopt_val <- log(x$opt_penalty)
      }
    }
    graphics::plot(
      x = xval,
      y = cverr,
      ylab = paste0("Mean CV Error (", x$loss, ")"),
      xlab = xlab,
      ylim = range(c(cverr - cvsd, cverr + cvsd)),
      type = "n"
    )
    graphics::arrows(
      xval,
      cverr - cvsd,
      xval,
      cverr + cvsd,
      length = 0.025,
      angle = 90,
      code = 3,
      col = "lightgray"
    )
    graphics::points(
      x = xval,
      y = cverr,
      col = "dodgerblue4",
      pch = 16,
    )
    graphics::abline(v = xopt_val, col = "firebrick")
  } else {
    cvgrid <- x$cv_mean
    cvgrid <- cvgrid[rev(seq_len(nrow(cvgrid))), ]
    cvgrid <- cvgrid[, rev(seq_len(ncol(cvgrid)))]
    minx <- log(x$opt_penalty_ext)
    miny <- log(x$opt_penalty)

    contour_colors <- c(
      "#014636", "#016C59", "#02818A", "#3690C0",
      "#67A9CF", "#A6BDDB", "#D0D1E6", "#ECE2F0", "#FFF7FB"
    )

    graphics::filled.contour(
      x = log(as.numeric(colnames(cvgrid))),
      y = log(as.numeric(rownames(cvgrid))),
      z = t(cvgrid),
      col = colorRampPalette(contour_colors)(25),
      xlab = "log(External Penalty)",
      ylab = "log(Penalty)",
      plot.axes = {
        axis(1)
        axis(2)
        points(minx, miny, col = "red", pch = 16)
      }
    )
  }
}
