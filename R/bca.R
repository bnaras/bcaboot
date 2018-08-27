## Version May 30, 2018

.cou <- function (x, b, labels = 0) {
    b0 <- min(x) - 1e-05
    b1 <- max(x) + 1e-05
    if (b0 < min(b))
        b <- c(b0, b)
    if (b1 > max(b))
        b <- c(b, b1)
    v <- table(cut(x, b))
    if (labels == 0) {
        as.vector(v)
    } else {
        v
    }
}

bca <- function(t0, tt, bb, ASL,
                alpha = c(0.025, 0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95, 0.975),
                trun = 0.001, pct = 0.5, ttrun = 0.005, zc = 7.5, sw = 0) {
    ## this version of bca accepts trun=vector; also returns amat

    # if bb Bxp matrix of bootsims suff vector corresponding to tt then a is computed
    # from skewness of local stats 'd'

    call <- sys.call()

    if (!missing(bb)) {
        bb <- as.matrix(bb)
        W <- scale(bb)
        S <- as.vector((W^2) %*% rep(1, ncol(W)))
        s <- quantile(S, pct)
        ii <- 1:length(tt)
        ii <- ii[S < s]
        t. <- stats::lm(tt[ii] ~ bb[ii, ])$coef[-1]
        bb <- scale(bb, T, F)
        s. <- stats::lm(tt ~ bb)$coef[-1]
        u. <- 2 * t. - s.
        V <- stats::cov(bb)
        sdu <- drop((t(u.) %*% V %*% u.)^0.5)
        if (sw == 2)
            return(list(t. = t., ii = ii, sdu = sdu))
        d <- as.vector(bb %*% t.)
        ntrun <- length(trun)

        amat <- matrix(0, ntrun, 3)
        dimnames(amat) <- list(trun, c("a", "az", "sdd"))
        for (k in 1:ntrun) {
            trunk <- trun[k]
            ulo <- quantile(d, trunk)
            uup <- quantile(d, 1 - trunk)
            dk <- pmax(pmin(d, uup), ulo)
            md <- mean(dk)
            sdk <- stats::sd(dk)
            ak <- (1/6) * mean((dk - md)^3)/sdk^3
            dz <- 1/(1 + exp(pmin(pmax(zc * d/sdk, -100), 100)))

            akz <- stats::qnorm(mean(dz))
            amat[k, ] <- c(ak, akz, sdk)
        }
    }
    trun <- trun[1]
    sdd <- amat[1, 3]
    a <- amat[1, 1]
    az <- amat[1, 2]
    if (sw == 3)
        return(d)

    Stand <- t0 + stats::sd(tt) * stats::qnorm(alpha)
    B <- length(tt)
    sd.mean.B <- c(stats::sd(tt), sdd, mean(tt), B)

    z0 <- stats::qnorm(mean(tt < t0))
    Z <- z0 + stats::qnorm(alpha)

    atil <- stats::pnorm(z0 + Z/(1 - a * Z))
    lims <- quantile(tt, atil)
    lims <- rbind(lims, atil, Stand)
    dimnames(lims) <- list(c("bcalims", "%iles", "Stand"), alpha)
    lims <- t(lims)
    if (missing(t0)) {
        thetahat <- NA
    } else {
        thetahat <- t0
    }

    A <- NA
    if (!missing(bb)) {
        A <- mean((tt - mean(tt))^2 * d)/(2 * var(d) * stats::sd(tt))
    }
    thet.a.z0.A <- c(thetahat, a, z0, A, az)
    u0 <- 2 * t0 - mean(tt)
    ustats <- c(u0, sdu)
    vl <- list(call, lims, thet.a.z0.A, sd.mean.B, ustats)
    names(vl) <- c("call", "lims", "thet.a.z0.A.az", "sd.sdd.mean.B", "ustats")
    if (!missing(ASL)) {
        z1 <- cumsum(.cou(tt, ASL))/length(tt)
        z1 <- stats::qnorm(z1[-(length(z1))]) - z0
        asl <- stats::pnorm(z1/(1 + a * z1) - z0)
        asl <- rbind(ASL, asl)
        vl$asl <- asl
    }
    if (ntrun > 1)
        vl$amat <- amat
    vl
}
