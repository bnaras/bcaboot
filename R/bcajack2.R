## Version of June 16, 2018

#' @title Nonparametric bias-corrected and accelerated bootstrap
#'     confidence limits

#' @description This function is a version of `bcajack` that allows
#'     all the recomputations of the original statistic function
#'     \eqn{f} to be carried out separately. This is an advantage
#'     if \eqn{f} is time-consuming, in which case the B
#'     replications for the nonparametric bca calculations might need
#'     to be done on a distributed basis.
#'
#' To use `bcajack2` in this mode, we first compute a list `Blist` via
#' `Blist <- list(Y = Y, tt = tt, t0 = t0)`.  Here `tt` is a vector of
#' length `B` having i-th entry `tt[i] <- func(x[Ii,], ...)`, where `x`
#' is the \eqn{n \times p} data matrix and `Ii` is a bootstrap vector
#' of (observation) indices. `Y` is a `B` by \eqn{n} count matrix,
#' whose i-th row is the counts corresponding to `Ii`. For example if
#' n = 5 and `Ii = (2, 5, 2, 1, 4)`, then `Yi = (1, 2, 0, 1,
#' 1)`. Having computed `Blist`, `bcajack2` is invoked as
#' `bcajack2(Blist)` without need to enter the function \eqn{func}.
#'
#' @inheritParams bcajack
#'
#' @param B number of bootstrap replications. `B` can also be a vector
#'     of `B` bootstrap replications of the estimated parameter of
#'     interest, computed separately. If `B` is `Blist` as explained
#'     above, `x` is not needed.
#' @param pct `bcajack2` uses those count vectors nearest (1,1,...1)
#'     to estimate the gradient of the statistic, "nearest" being
#'     defined as those count vectors in the smallest `pct` of all B
#'     of them. Default value for `pct is 1/3 (see appendix in Efron
#'     and Narasimhan for further details)
#'
#' @return a named list of several items
#'
#' * __lims__ : first column shows the estimated bca confidence limits
#'     at the requested alpha percentiles. These can be compared with
#'     the standard limits \eqn{\hat{\theta} +
#'     \hat{\sigma}z_{\alpha}}, third column. The second column
#'     `jacksd` gives the internal standard errors for the bca limits,
#'     quite small in the example. Column 4, `pct`, gives the
#'     percentiles of the ordered B bootstrap replications
#'     corresponding to the bca limits, eg the 897th largest
#'     replication equalling the .975 bca limit .557.
#'
#' * __stats__ : top line of stats shows 5 estimates: theta is
#'     \eqn{func(x)}, original point estimate of the parameter of
#'     interest; `sdboot` is its bootstrap estimate of standard error;
#'     `z0` is the bca bias correction value, in this case quite
#'     negative; `a` is the _acceleration_, a component of the bca
#'     limits (nearly zero here); `sdjack` is the jackknife estimate
#'     of standard error for theta. Bottom line gives the internal
#'     standard errors for the five quantities above. This is
#'     substantial for `z0` above.
#'
#' * __B.mean__ : bootstrap sample size B, and the mean of the B
#'     bootstrap replications \eqn{\hat{\theta^*}}
#'
#' * __ustats__ : The bias-corrected estimator `2 * t0 - mean(tt)`,
#'     and an estimate `sdu` of its sampling error
#'
#' * __seed__ : The random number state for reproducibility
#'
#' @importFrom stats quantile
#' @export
#' @examples
#' data(diabetes, package = "bcaboot")
#' Xy <- cbind(diabetes$x, diabetes$y)
#' rfun <- function(Xy) {
#'   y <- Xy[, 11]
#'   X <- Xy[, 1:10]
#'   summary(lm(y~X) )$adj.r.squared
#' }
#' set.seed(1234)
#' bcajack2(x = Xy, B = 1000, func = rfun, m = 40, verbose = FALSE)
#'
#' @export
bcajack2 <- function(x, B, func, ..., m = nrow(x), mr, pct = 0.333, K = 2, J = 12,
                     alpha = c(0.025, 0.05, 0.1, 0.16),
                     verbose = TRUE) {

    call <- match.call()

    ## Save rng state
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        stats::runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    qbca2 <- function(Y, tt, t0, alpha, pct) {
        m <- ncol(Y)
        B <- nrow(Y)
        o1 <- rep(1, m)
        D <- rep(0, B)

        for (i in seq_len(B)) {
            Yi <- Y[i, ]
            d <- 2 * Yi * log(Yi) - 2 * (Yi - 1)
            d[Yi == 0] <- 2
            D[i] <- sum(d)
        }
        Qd <- stats::quantile(D, pct)
        ip <- seq_len(B)[D <= Qd]
        ty. <- as.vector(m * stats::lm(tt[ip] ~ Y[ip, ] - 1)$coef)
        ty. <- ty. - mean(ty.)
        a <- (1/6) * sum(ty.^3)/sum(ty.^2)^1.5
        ## if (sw == 3)
        ##     return(ty.)
        s <- mean(tt)
        B.mean <- c(B, s)

        zalpha <- stats::qnorm(alpha)
        nal <- length(alpha)
        ustat <- 2 * t0 - s
        ##s. <- m * .v(stats::cov(tt, Y))
        s. <- m * as.vector(stats::cov(tt, Y))
        u. <- 2 * ty. - s.
        sdu <- sum(u.^2)^0.5/m
        ustats <- c(ustat, sdu)
        names(ustats) <- c("ustat", "sdu")

        sdboot <- stats::sd(tt)
        sdjack <- sqrt(sum(ty.^2))/(m - 1)
        z0 <- stats::qnorm(sum(tt < t0)/B)

        iles <- stats::pnorm(z0 + (z0 + zalpha)/(1 - a * (z0 + zalpha)))
        ooo <- trunc(iles * B)
        ooo <- pmin(pmax(ooo, 1), B)
        lims <- sort(tt)[ooo]
        standard <- t0 + sdboot * stats::qnorm(alpha)
        ##lims <- round(cbind(lims, standard), rou)
        lims <- cbind(lims, standard)
        dimnames(lims) <- list(alpha, c("bca", "std"))
        stats <- c(t0, sdboot, z0, a, sdjack)
        names(stats) <- c("theta", "sdboot", "z0", "a", "sdjack")
        vl <- list(lims = lims, stats = stats, B.mean = B.mean, ustats = ustats)
        return(vl)
    }

    alpha <- alpha[alpha < 0.5]
    alpha <- c(alpha, 0.5, rev(1 - alpha))

    if (is.list(B)) {
        Y <- B$Y
        tt <- B$tt
        t0 <- B$t0
        B <- length(tt)
        ##vl0 <- qbca2(Y, tt, t0, alpha = alpha, pct = pct, rou = rou)
        vl0 <- qbca2(Y, tt, t0, alpha = alpha, pct = pct)
    } else {
        if (is.vector(x))
            x <- as.matrix(x)
        n <- nrow(x)
        tt <- rep(0, B)
        t0 <- func(x, ...)

        if (m == n) {
            ii <- sample(x = seq_len(n), size = n * B, replace = TRUE)
            ii <- matrix(ii, B)
            Y <- matrix(0, B, n)
            if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
            for (k in 1:B) {
                ik <- ii[k, ]
                tt[k] <- func(x[ik, ], ...)
                Y[k, ] <- table(c(ik, 1:n)) - 1
                if (verbose) utils::setTxtProgressBar(pb, k)
            }
            if (verbose) close(pb)
            ##vl0 <- qbca2(Y, tt, t0, alpha = alpha, pct = pct, rou = rou)
            vl0 <- qbca2(Y, tt, t0, alpha = alpha, pct = pct)
        }

        if (m < n) {
            r <- n%%m
            Imat <- matrix(sample(x = seq_len(n), size = n - r), m)
            Iout <- setdiff(1:n, Imat)
            ii <- sample(x = seq_len(m), size = m * B, replace = TRUE)
            ii <- matrix(ii, B)
            Y <- matrix(0, B, m)
            if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
            for (k in 1:B) {
                ik <- ii[k, ]
                Ik <- c(t(Imat[ik, ]))
                Ik <- c(Ik, Iout)
                tt[k] <- func(x[Ik, ], ...)
                Y[k, ] <- table(c(ik, 1:m)) - 1
                if (verbose) utils::setTxtProgressBar(pb, k)
            }
            if (verbose) close(pb)
            ##vl0 <- qbca2(Y, tt, t0, alpha = alpha, pct = pct, rou = rou)
            vl0 <- qbca2(Y, tt, t0, alpha = alpha, pct = pct)
        }
    }
    vl0$call <- call
    vl0$seed <- seed

    if (K == 0)
        bcaboot.return(vl0)

    nal <- length(alpha)
    Pct <- rep(0, nal)
    ##for (i in 1:nal) Pct[i] <- round(sum(tt <= vl0$lims[i, 1])/B, rou)
    for (i in 1:nal) Pct[i] <- sum(tt <= vl0$lims[i, 1])/B
    Stand <- vl0$stats[1] + vl0$stats[2] * stats::qnorm(alpha)
    Limsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)

    Limbcsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)
    for (k in 1:K) {
        rr <- B%%J
        II <- sample(x = B, size = B - rr)
        II <- matrix(II, ncol = J)
        limbc <- limst <- matrix(0, length(alpha), J)
        stats <- matrix(0, 5, J)

        for (j in 1:J) {
            iij <- c(II[, -j])
            Yj <- Y[iij, ]
            ttj <- tt[iij]
            ##vlj <- qbca2(Yj, ttj, t0, alpha, pct, rou)
            vlj <- qbca2(Yj, ttj, t0, alpha, pct)
            limbc[, j] <- vlj$lims[, 1]
            limst[, j] <- vlj$lims[, 2]
            stats[, j] <- vlj$stats
        }

        ## if (sw == 4)
        ##     return(list(limbc = limbc, limst = limst, stats = stats))
        Limbcsd[, k] <- apply(limbc, 1, sd) * (J - 1)/sqrt(J)
        Statsd[, k] <- apply(stats, 1, sd) * (J - 1)/sqrt(J)
        ## if (verbose)
        ##     cat("{", k, "}", sep = "")
        ## if (sw == 6)
        ##     return(list(Limbcsd = Limbcsd, Statsd = Statsd))
    }
    limsd <- rowMeans(Limbcsd, 1)
    statsd <- rowMeans(Statsd, 1)
    ##limits <- round(cbind(vl0$lims[, 1], limsd, vl0$lims[, 2], Pct), rou)
    limits <- cbind(vl0$lims[, 1], limsd, vl0$lims[, 2], Pct)
    dimnames(limits) <- list(alpha, c("bca", "jacksd", "std", "pct"))
    ##stats <- round(rbind(vl0$stats, statsd), rou)
    stats <- rbind(vl0$stats, statsd)
    ##ustats <- round(vl0$ustats, rou)
    ustats <- vl0$ustats
    ##B.mean <- c(B, round(mean(tt), rou))
    B.mean <- c(B, mean(tt))
    dimnames(stats) <- list(c("est", "jsd"), c("theta", "sdboot", "z0", "a", "sdjack"))
    vll <- list(call = call, lims = limits, stats = stats, B.mean = B.mean, ustats = ustats,
                seed = seed)
    ## if (sw == 5)
    ##     vll$tt <- tt
    bcaboot.return(vll)
}


