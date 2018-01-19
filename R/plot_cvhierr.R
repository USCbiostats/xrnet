#' Plot k-fold cross-validation error grid
#'
#' @param x A cvhierr class object
#' @param ... Additional graphics parameters

#' @export
#' @importFrom graphics filled.contour axis points
#' @importFrom grDevices colorRampPalette

plot.cvhierr <- function(x, ...) {
    cvgrid <- x$cv_mean
    cvgrid <- cvgrid[rev(seq_len(nrow(cvgrid))), ]
    cvgrid <- cvgrid[ , rev(seq_len(ncol(cvgrid)))]
    minx <- log(x$minl2)
    miny <- log(x$minl1)

    contour_colors <- c("#014636", "#016C59", "#02818A", "#3690C0",
                        "#67A9CF", "#A6BDDB", "#D0D1E6", "#ECE2F0", "#FFF7FB")

    filled.contour(x = log(as.numeric(colnames(cvgrid))),
                   y = log(as.numeric(rownames(cvgrid))),
                   z = t(cvgrid),
                   col = colorRampPalette(contour_colors)(25),
                   xlab = expression(log(lambda[2])),
                   ylab = expression(log(lambda[1])),
                   plot.axes = {axis(1); axis(2); points(minx, miny, col = "red", pch = 16)})

}
