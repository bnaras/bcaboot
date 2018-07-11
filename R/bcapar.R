## Version of May 30, 2018

bcapar <- function(t0, tt, bb, al = c(0.025, 0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95,
    0.975), J = 10, K = 6, trun = 0.001, pct = 0.333, roun = 3, cd = 0, func) {
    # if(!missing(func)) adds abc limits and stats K=0 skips jackknifing computes
    # confidence density weights 'w' if cd=1

    if (K == 0)
        return(bca(t0, tt, bb, al = al, trun = trun, pct = pct))
    syscall <- match.call()
    B <- length(tt)
    vl0 <- bca(t0, tt, bb = as.matrix(bb), al = al, trun = trun, pct = pct)
    Stand <- vl0$thet[1] + vl0$sd.m[1] * qnorm(al)
    Limsd <- matrix(0, length(al), K)
    Thsd <- matrix(0, 5, K)
    Sdmsd <- matrix(0, 4, K)
    Ustm <- matrix(0, 2, K)
    for (k in 1:K) {
        II <- sample(B, B)
        II <- matrix(II, ncol = J)
        lims <- matrix(0, length(al), J)
        th <- matrix(0, 5, J)
        sdm <- matrix(0, 4, J)
        usm <- matrix(0, 2, J)
        for (j in 1:J) {
            iij <- c(II[, -j])
            bbj <- as.matrix(bb[iij, ])
            ttj <- tt[iij]
            vlj <- bca(t0, ttj, bb = bbj, al = al, trun = trun, pct = pct)
            lims[, j] <- vlj$lims[, 1]
            th[, j] <- vlj$thet.a
            sdm[, j] <- vlj$sd.s
            usm[, j] <- vlj$ustats
        }

        Limsd[, k] <- apply(lims, 1, sd) * (J - 1)/sqrt(J)
        Thsd[, k] <- apply(th, 1, sd) * (J - 1)/sqrt(J)
        Sdmsd[, k] <- apply(sdm, 1, sd) * (J - 1)/sqrt(J)
        Ustm[, k] <- apply(usm, 1, sd) * (J - 1)/sqrt(J)
    }

    limsd <- rowMeans(Limsd)
    thsd <- rowMeans(Thsd)
    sdmsd <- rowMeans(Sdmsd)
    ustsd <- rowMeans(Ustm)
    lim0 <- vl0$lims
    lim0 <- cbind(lim0[, 1], limsd, lim0[, 2:3])
    dimnames(lim0) <- list(al, c("bcalims", "jacksd", "pctiles", "Stand"))
    th0 <- vl0$thet
    th0 <- rbind(th0, thsd)
    dimnames(th0) <- list(c("est", "jsd"), c("thet", "a", "z0", "A", "az"))
    sdm0 <- vl0$sd.s
    sdm0 <- rbind(sdm0, sdmsd)
    dimnames(sdm0) <- list(c("est", "jsd"), c("sd", "sdd", "mean", "B"))
    ust0 <- vl0$ustats
    ust0 <- rbind(ust0, ustsd)
    ust0 <- cbind(ust0, c(B, 0))
    dimnames(ust0) <- list(c("est", "jsd"), c("ustat", "sdu", "B"))
    lims <- round(lim0, roun)
    th0 <- round(th0, roun)
    sdm0 <- round(sdm0, roun)
    ust0 <- round(ust0, roun)
    stats <- cbind(th0, sdm0)
    stats <- stats[, c(1, 6, 2, 5, 3, 4, 7, 8)]

    vl <- list(syscall = syscall, lims = lims, stats = stats, ustats = ust0)
    if (length(trun) > 1)
        vl$amat <- vl0$amat
    if (cd == 1) {
        a <- stats[1, 4]
        z0 <- stats[1, 5]
        G <- (rank(tt) - 0.5)/B
        zth <- qnorm(G) - z0
        az <- 1 + a * zth
        num <- dnorm(zth/az - z0)
        den <- az^2 * dnorm(zth + z0)
        w <- num/den
        vl$w <- w
    }

    if (!missing(func)) {
        vla <- abcpar(func, bb, alpha = al[al < 0.5])
        abclims <- vla$lims[, 1]
        vl$lims <- round(cbind(vl$lims, abclims), roun)
        abcstats <- round(vla$a.z0.cq[1:2], roun)
        names(abcstats) <- c("a", "z0")
        vl$abcstats <- abcstats
    }
    vl
}
