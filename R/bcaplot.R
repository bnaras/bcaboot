#' @title Plots of bca confidence limits
#'
#' @description \code{bcaplot} uses the output of \code{bcajack},
#'     \code{bcajack2}, or \code{bcapar} to plot bca and standard
#'     confidence limits for the parameter of interest.
#'
#' @param vl output of \code{bcajack}, \code{bcajack2}, or
#'     \code{bcapar}
#' @param main The main caption (can be empty)
#' @param xlab The x axis label (supplied if not specified)
#' @param ylab The y axis labels (supplied if not specified)
#' @param alpha Coverages are \eqn{1-2\alpha},
#'     e.g. \code{alpha=c(.025,.05)} plots intervals
#'     \code{[.025,.975]} and \code{[.05,.95]}. Default is
#'     \code{alpha=c(.025,.05,.1,.16)} giving coverages
#'     .95,.90,.80,.68
#' @param ylim y axis plot limits set automatically if not provided
#' @param xlim x axis plot limits set automatically if not provided
#' @param add \code{add=1} adds a new plot of bca limits (in red) to
#'     an existing plot
#' @param sub subtitle (can be empty)
#' @param sw \code{sw=1} draws light vertical dashed lines showing the
#'     bca intervals
#' @param ... further args for plot
#'
#' @details confidence interval endpoints are plotted vertically versus
#'     two-sided coverages \eqn{1-2\alpha}. Bca limits in black,
#'     Standard limits in green (dashed.). If \code{vl$lims} includes
#'     the column "jacksd" of jackknife internal standard deviations
#'     then these are indicated by vertical red bars centered at the
#'     bca limit points.
#' @export
#'
bcaplot <- function(vl, main = " ", xlab = "coverage", ylab = "limits", alpha = c(0.025,
    0.05, 0.1, 0.16), ylim, xlim, add = 0, sub = "black=bca, green=standard", sw = 1,
    ...) {
    alp <- rev(1 - 2 * alpha)
    vm <- vl$lims

    if (ncol(vm) == 2) {
        li <- vl$lims
        o <- rep(0, nrow(li))
        li <- cbind(li[, 1], o, li[, 2])
        dimnames(li)[[2]] <- c("bca", "jacksd", "std")
        sta <- vl$stats
        sta <- rbind(sta, 0)
        vl$stats <- sta
    }
    if (nrow(vm) == 9)
        vm <- vm[-5, ]

    vn <- dimnames(vm)[[2]]
    nvn <- length(vn)
    jn <- (1:nvn)[vn == "std"]
    vja <- (1:nvn)[vn == "jacksd"]
    if (length(vja) > 0)
        jasd <- vm[, vja]

    vm <- vm[, c(1, jn)]
    if (add == 1) {
        graphics::lines(alp, vm[4:1, 1], lwd = 2, lty = 2, col = 2)
        graphics::lines(alp, vm[5:8, 1], lwd = 2, lty = 2, col = 2)
    }
    if (add == 0) {
        if (length(vl$thet) == 0)
            thet <- vl$stat[1, 1] else thet <- vl$thet[1]
        if (missing(ylim))
            ylim <- range(vm)
        if (missing(xlim))
            xlim <- c(min(alp), max(alp) + 0.005)
        graphics::plot(alp, vm[4:1, 1], type = "l", lwd = 3, ..., ylim = ylim, xlab = xlab,
                       ylab = ylab, main = main, xlim = xlim, sub = sub)
        graphics::lines(alp, vm[5:8, 1], lwd = 3)
        graphics::lines(alp, vm[4:1, 2], lwd = 2, lty = 3, col = 3)
        graphics::lines(alp, vm[5:8, 2], lwd = 3, lty = 3, col = 3)
        graphics::abline(thet, 0, col = 2, lty = 2)
        graphics::points(alp, vm[4:1, 1], pch = 16, cex = 1.1)
        graphics::points(alp, vm[5:8, 1], pch = 16, cex = 1.1)
        for (k in 1:length(alp)) graphics::text(alp[k], thet, paste(100 * alp[k], "%", sep = ""),
            cex = 1.2)
        if (length(vja) > 0) {
            cover <- c(0.95, 0.9, 0.8, 0.68, 0.68, 0.8, 0.9, 0.95)
            graphics::segments(cover, vm[, 1] - jasd, cover, vm[, 1] + jasd, col = 2, lty = 1,
                               lwd = 3)
        }
        if (sw == 1)
            graphics::segments(alp, vm[4:1, 1], alp, vm[5:8, 1], lwd = 1, lty = 2, col = 1)
    }
}
