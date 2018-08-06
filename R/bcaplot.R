#' @export
bcaplot <- function(vl, main = " ", xlab = "coverage", ylab = "limits", alpha = c(0.025,
    0.05, 0.1, 0.16), ylim, xlim, add = 0, sub = "black=bca, green=standard", sw = 1,
    ...) {
    alp <- rev(1 - 2 * alpha)
    vm <- vl$lims

    if (ncol(vm) == 2) {
        li <- vl$lims
        o <- rep(0, nrow(li))
        li <- cbind(li[, 1], o, li[, 2])
        dimnames(li)[[2]] <- c("bcalims", "jacksd", "standard")
        sta <- vl$stats
        sta <- rbind(sta, 0)
        vl$stats <- sta
    }
    if (nrow(vm) == 9)
        vm <- vm[-5, ]

    vn <- dimnames(vm)[[2]]
    nvn <- length(vn)
    jn <- (1:nvn)[vn == "standard" | vn == "Stand"]
    vja <- (1:nvn)[vn == "jacksd"]
    if (length(vja) > 0)
        jasd <- vm[, vja]

    vm <- vm[, c(1, jn)]
    if (add == 1) {
        lines(alp, vm[4:1, 1], lwd = 2, lty = 2, col = 2)
        lines(alp, vm[5:8, 1], lwd = 2, lty = 2, col = 2)
    }
    if (add == 0) {
        if (length(vl$thet) == 0)
            thet <- vl$stat[1, 1] else thet <- vl$thet[1]
        if (missing(ylim))
            ylim <- range(vm)
        if (missing(xlim))
            xlim <- c(min(alp), max(alp) + 0.005)
        plot(alp, vm[4:1, 1], type = "l", lwd = 3, ..., ylim = ylim, xlab = xlab,
            ylab = ylab, main = main, xlim = xlim, sub = sub)
        lines(alp, vm[5:8, 1], lwd = 3)
        lines(alp, vm[4:1, 2], lwd = 2, lty = 3, col = 3)
        lines(alp, vm[5:8, 2], lwd = 3, lty = 3, col = 3)
        abline(thet, 0, col = 2, lty = 2)
        points(alp, vm[4:1, 1], pch = 16, cex = 1.1)
        points(alp, vm[5:8, 1], pch = 16, cex = 1.1)
        for (k in 1:length(alp)) text(alp[k], thet, paste(100 * alp[k], "%", sep = ""),
            cex = 1.2)
        if (length(vja) > 0) {
            cover <- c(0.95, 0.9, 0.8, 0.68, 0.68, 0.8, 0.9, 0.95)
            segments(cover, vm[, 1] - jasd, cover, vm[, 1] + jasd, col = 2, lty = 1,
                lwd = 3)
        }
        if (sw == 1)
            segments(alp, vm[4:1, 1], alp, vm[5:8, 1], lwd = 1, lty = 2, col = 1)
    }
}
