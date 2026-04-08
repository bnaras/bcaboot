## Version May 30, 2018
##
## Internal parametric BCa engine. Called by bcapar().
## Computes acceleration from skewness of projected direction d = bb %*% local_grad,
## where local_grad is the gradient of the statistic estimated via local regression
## on the sufficient statistics bb.

## Bin-counting utility for ASL (achieved significance level) computation.
bin_counts <- function (x, b, labels = 0) {
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
                alpha, trun = 0.001, pct = 0.5, ttrun = 0.005, zc = 7.5, sw = 0) {

    call <- sys.call()
    stopifnot(!missing(bb))
    bb <- as.matrix(bb)
    ## S = Mahalanobis-like distance: sum of squared standardized components.
    ## Select the pct fraction of sufficient statistic vectors nearest to center.
    W <- scale(bb)
    S <- as.vector((W^2) %*% rep(1, ncol(W)))
    s <- quantile(S, pct)
    ii <- seq_along(tt)
    ii <- ii[S < s]
    ## local_grad = gradient of statistic via regression on nearby bb vectors
    local_grad <- stats::lm(tt[ii] ~ bb[ii, ])$coef[-1]
    bb <- scale(bb, center = TRUE, scale = FALSE)
    ## global_grad = gradient via regression on all bb vectors
    global_grad <- stats::lm(tt ~ bb)$coef[-1]
    ## Variance-stabilized gradient: 2*local - global. EN20 Sec 3
    adj_grad <- 2 * local_grad - global_grad
    V <- stats::cov(bb)
    sdu <- drop((t(adj_grad) %*% V %*% adj_grad)^0.5)
    if (sw == 2)
        return(list(local_grad = local_grad, ii = ii, sdu = sdu))
    ## d = projection of centered bb onto gradient direction
    d <- as.vector(bb %*% local_grad)
    ntrun <- length(trun)

    ## Compute acceleration at each truncation level.
    ## Truncation clips extreme values of d to reduce outlier sensitivity.
    ## ak = standard acceleration via skewness of truncated d. Efron (1987) Sec 6.
    ## akz = alternative acceleration via logistic transform (more robust). EN20 Sec 4.
    ## sdd = delta-method SE estimate from projection.
    amat <- matrix(0, ntrun, 3)
    dimnames(amat) <- list(trun, c("a", "az", "sdd"))
    for (k in seq_len(ntrun)) {
        trunk <- trun[k]
        ulo <- quantile(d, trunk)
        uup <- quantile(d, 1 - trunk)
        dk <- pmax(pmin(d, uup), ulo)
        md <- mean(dk)
        sdk <- stats::sd(dk)
        ak <- (1/6) * mean((dk - md)^3)/sdk^3
        ## Sigmoid-based acceleration: maps d to (0,1) via logistic,
        ## then inverts through Phi. More robust to extreme d. EN20 Sec 4.
        dz <- 1/(1 + exp(pmin(pmax(zc * d/sdk, -100), 100)))

        akz <- stats::qnorm(mean(dz))
        amat[k, ] <- c(ak, akz, sdk)
    }
    trun <- trun[1]
    sdd <- amat[1, 3]
    a <- amat[1, 1]
    az <- amat[1, 2]
    if (sw == 3)
        return(d)

    Stand <- t0 + stats::sd(tt) * stats::qnorm(alpha)
    B <- length(tt)

    ## Bias-correction z0. Efron (1987) Sec 2
    z0 <- stats::qnorm(mean(tt < t0))
    zalpha <- stats::qnorm(alpha)

    ## BCa percentile formula. Efron (1987) Sec 2
    bca_result <- compute_bca_limits(z0, a, zalpha, tt)
    atil <- bca_result$iles
    lims <- quantile(tt, atil)
    lims <- rbind(lims, atil, Stand)
    dimnames(lims) <- list(c("bcalims", "%iles", "Stand"), alpha)
    lims <- t(lims)
    thetahat <- if (missing(t0)) NA else t0

    ## Big-A raw acceleration: E[(tt-mean)^2 * d] / (2*Var(d)*sd(tt)).
    ## Delta-method correction measuring rate of SE change. Efron (1987) Sec 7.
    A <- mean((tt - mean(tt))^2 * d)/(2 * var(d) * stats::sd(tt))
    theta_a_z0_A_az <- c(thetahat, a, z0, A, az)
    names(theta_a_z0_A_az) <- c("theta", "a", "z0", "A", "az")
    ## Bias-corrected point estimate: 2*t0 - mean(tt)
    u0 <- 2 * t0 - mean(tt)
    ustats <- c(u0, sdu)
    sd_sdd_mean_B <- c(stats::sd(tt), sdd, mean(tt), B)
    names(sd_sdd_mean_B) <- c("sd", "sdd", "mean", "B")
    result <- list(call = call, lims = lims,
                   theta_a_z0_A_az = theta_a_z0_A_az,
                   sd_sdd_mean_B = sd_sdd_mean_B,
                   ustats = ustats)
    if (!missing(ASL)) {
        z1 <- cumsum(bin_counts(tt, ASL))/length(tt)
        z1 <- stats::qnorm(z1[-(length(z1))]) - z0
        asl <- stats::pnorm(z1/(1 + a * z1) - z0)
        asl <- rbind(ASL, asl)
        result$asl <- asl
    }
    if (ntrun > 1)
        result$amat <- amat
    result
}
