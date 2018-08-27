## Version of June 1, 2018

#'
#' @title Nonparametric bias-corrected and accelerated bootstrap
#'     confidence limits
#'
#' @description This routine computes nonparametric confidence
#'     intervals for bootstrap estimates. For reproducibility, save or
#'     set the random number state before calling this routine.
#'
#' @details
#' Bootstrap confidence intervals depend on three elements:
#'
#' - the cdf of the \eqn{B} bootstrap replications \eqn{t_i^*}, \eqn{i=1\ldots B}
#' - the bias-correction number \eqn{z_0=\Phi(\sum_i^B I(t_i^* < t_0) / B )}
#'   where \eqn{t_0=f(x)} is the original estimate
#' - the acceleration number \eqn{a} that measures the rate of
#'   change in \eqn{\sigma_{t_0}} as \eqn{x}, the data changes.
#'
#' The first two of these depend only on the bootstrap distribution,
#' and not how it is generated: parametrically or
#' non-parametrically. Program bcajack can be used in a hybrid fashion
#' in which the vector \code{tt} of B bootstrap replications is first
#' generated from a parametric model.
#'
#' So, in the diabetes example below, we might first draw bootstrap
#' samples \eqn{y^* \sim N(X\hat{\beta}, \hat{\sigma}^2 I)} where
#' \eqn{\hat{\beta}} and \eqn{\hat{\sigma}} were obtained from
#' `lm(y~X)`; each \eqn{y^*} would then provide a bootstrap
#' replication `tstar = rfun(cbind(X, ystar))`.  Then we could get bca
#' intervals from `bcajack(Xy, tt, rfun ....)` with `tt`,
#' the vector of B `tstar` values. The only difference from a full
#' parametric bca analysis would lie in the nonparametric estimation
#' of \eqn{a}, often a negligible error.
#'
#' @param x an \eqn{n \times p} data matrix, rows are observed
#'     \eqn{p}-vectors, assumed to be independently sampled from
#'     target population. If \eqn{p} is 1 then `x` can be a vector.
#' @param B number of bootstrap replications. It can also be a vector
#'     of `B` bootstrap replications of the estimated parameter of
#'     interest, computed separately.
#' @param func function \eqn{\hat{\theta}=func(x)} computing estimate of the
#'     parameter of interest; \eqn{func(x)} should return a real value
#'     for any \eqn{n^\prime \times p} matrix \eqn{x^\prime},
#'     \eqn{n^\prime} not necessarily equal to \eqn{n}
#' @param ... additional arguments for `func`.
#' @param m an integer less than or equal to \eqn{n}; the routine
#'     collects the \eqn{n} rows of `x` into `m` groups to speed up
#'     the jackknife calculations for estimating the acceleration
#'     value \eqn{a}; typically `m` is 20 or 40 and does not have to
#'     exactly divide \eqn{n}. However, warnings will be shown.
#' @param mr if \eqn{m < n} then `mr` repetions of the randomly
#'     grouped jackknife calculations are averaged.
#' @param K a non-negative integer. If `K` > 0, bcajack also returns
#'     estimates of _internal standard error_, that is, of the
#'     variability due to stopping at `B` bootstrap replications
#'     rather than going on to infinity. These are obtained from a
#'     second type of jackknifing, taking an average of `K` separate
#'     jackknife estimates, each randomly splitting the `B` bootstrap
#'     replications into `J` groups.
#' @param J the number of groups into which the bootstrap replications
#'     are split
#' @param alpha percentiles desired for the bca confidence limits. One
#'     only needs to provide `alpha` values below 0.5; the upper
#'     limits are automatically computed
#' @param verbose logical for verbose progress messages
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
#'     \eqn{f(x)}, original point estimate of the parameter of
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
#' @references DiCiccio T and Efron B (1996). Bootstrap confidence
#'     intervals. Statistical Science 11, 189-228
#' @references Efron B (1987). Better bootstrap confidence
#'     intervals. JASA 82 171-200
#' @references B. Efron and B. Narasimhan. Automatic Construction of
#'     Bootstrap Confidence Intervals, 2018.
#'
#' @importFrom graphics abline lines plot points segments text
#' @importFrom stats cov dnorm lm pnorm qnorm runif sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @examples
#' data(diabetes, package = "bcaboot")
#' Xy <- cbind(diabetes$x, diabetes$y)
#' rfun <- function(Xy) {
#'   y <- Xy[, 11]
#'   X <- Xy[, 1:10]
#'   summary(lm(y~X) )$adj.r.squared
#' }
#' set.seed(1234)
#' ## n = 442 = 34 * 13
#' bcajack(x = Xy, B = 1000, func = rfun, m = 34, verbose = FALSE)
#' @export
bcajack <- function(x, B, func, ..., m = nrow(x), mr = 5, K = 2, J = 10,
                    alpha = c(0.025, 0.05, 0.1, 0.16),  verbose = TRUE) {
    ## x is nxp data matrix, func is statistic thetahat=func(x) can enter #bootsize B
    ## for bootsim vector tt (which is calculated)

    call <- match.call()
    ## Save rng state
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        stats::runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    if (is.vector(x))
        x <- as.matrix(x)
    n <- nrow(x)
    nc <- ncol(x)
    ttind <- ifelse(length(B) > 1, 1, 0)
    if (ttind == 1) {
        tt <- B
        B <- length(tt)
    } else tt <- rep(0, B)
    t0 <- func(x, ...)
    u <- numeric(length = m)
    m1 <- sqrt(m * (m - 1))

    if (m == n) {
        for (i in seq_len(n)) {
            u[i] <- func(x[-i, ], ...)
        }
        t. <- (mean(u) - u) * (m - 1)
        a <- (1 / 6) * sum(t.^3) / (sum(t.^2))^1.5
        sdjack <- sqrt(sum(t.^2)) / m1
    }

    if (m < n) {
        ##aa <- ssj <- rep(0, mr)
        aa <- ssj <- numeric(mr)
        r <- n %% m
        seq_len_m <- seq_len(m)
        for (k in seq_len(mr)) {
            ##Imat <- matrix(sample(1:n, n - r), m)
            Imat <- sapply(seq_len_m, sample.int, n = n, size = n - r)
            Iout <- setdiff(seq_len(n), Imat)
            for (j in seq_len_m) {
                Ij <- setdiff(seq_len_m, j)
                ij <- c(c(Imat[Ij, ], Iout))
                u[j] <- func(x[ij, ])
            }
            t. <- (mean(u) - u) * (m - 1)
            aa[k] <- (1/6) * sum(t.^3)/(sum(t.^2))^1.5
            ssj[k] <- sqrt(sum(t.^2))/m1
        }
        a <- mean(aa)
        sdjack <- mean(ssj)
    }

    if (ttind == 0) {
        tY. <- Y. <- rep(0, n)
        if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
        for (j in seq_len(B)) {
            ij <- sample(x = n, size = n, replace = TRUE)
            Yj <- table(c(ij, 1:n)) - 1
            tt[j] <- func(x[ij, ], ...)
            tY. <- tY. + tt[j] * Yj
            Y. <- Y. + Yj
            if (verbose) utils::setTxtProgressBar(pb, j)
        }
        if (verbose) close(pb)
        tt. <- mean(tt)
        tY. <- tY./B
        Y. <- Y./B
        s. <- n * (tY. - tt. * Y.)
        u. <- 2 * t. - s.
        sdu <- sqrt(sum(u.^2))/n
        ustat <- 2 * t0 - tt.
        ustats <- c(ustat, sdu)
        names(ustats) <- c("ustat", "sdu")
    }
    B.mean <- c(B, mean(tt))
    alpha <- alpha[alpha < 0.5]
    alpha <- c(alpha, 0.5, rev(1 - alpha))

    zalpha <- stats::qnorm(alpha)
    nal <- length(alpha)

    sdboot0 <- stats::sd(tt)  # sdd=stats::sd(dd)
    z00 <- stats::qnorm(sum(tt < t0)/B)

    iles <- stats::pnorm(z00 + (z00 + zalpha)/(1 - a * (z00 + zalpha)))
    ooo <- trunc(iles * B)
    ooo <- pmin(pmax(ooo, 1), B)
    lims0 <- sort(tt)[ooo]
    standard <- t0 + sdboot0 * stats::qnorm(alpha)
    ## lims0 <- round(cbind(lims0, standard), rou)
    lims0 <- cbind(lims0, standard)
    dimnames(lims0) <- list(alpha, c("bca", "std"))
    ## stats0 <- round(c(t0, sdboot0, z00, a, sdjack), rou)
    stats0 <- c(t0, sdboot0, z00, a, sdjack)
    names(stats0) <- c("theta", "sdboot", "z0", "a", "sdjack")
    vl0 <- list(lims = lims0, stats = stats0, B.mean = B.mean, call = call, seed = seed)
    if (K == 0)
        bcaboot.return(vl0)

    pct <- rep(0, nal)
    ##for (i in 1:nal) pct[i] <- round(sum(tt <= lims0[i, 1])/B, 3)
    for (i in 1:nal) pct[i] <- sum(tt <= lims0[i, 1])/B
    Stand <- vl0$stats[1] + vl0$stats[2] * stats::qnorm(alpha)
    Limsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)

    for (k in 1:K) {
        II <- sample(x = B, size = B)
        II <- matrix(II, ncol = J)
        lims <- matrix(0, length(alpha), J)
        stats <- matrix(0, 5, J)
        for (j in 1:J) {
            iij <- c(II[, -j])
            ttj <- tt[iij]
            Bj <- length(ttj)
            sdboot <- stats::sd(ttj)
            z0 <- stats::qnorm(sum(ttj < t0)/Bj)

            iles <- stats::pnorm(z0 + (z0 + zalpha)/(1 - a * (z0 + zalpha)))
            oo <- trunc(iles * Bj)
            oo <- pmin(pmax(oo, 1), Bj)
            li <- sort(ttj)[oo]
            standard <- t0 + sdboot * stats::qnorm(alpha)
            ##sta <- round(c(t0, sdboot, z0, a, sdjack), rou)
            sta <- c(t0, sdboot, z0, a, sdjack)
            names(sta) <- c("theta", "sdboot", "z0", "a", "sdjack")
            lims[, j] <- li
            stats[, j] <- sta
        }
        Limsd[, k] <- apply(lims, 1, sd) * (J - 1)/sqrt(J)
        Statsd[, k] <- apply(stats, 1, sd) * (J - 1)/sqrt(J)
        ##if (verbose) cat("{", k, "}", sep = "")
    }
    limsd <- rowMeans(Limsd, 1)
    statsd <- rowMeans(Statsd, 1)
    ##limits <- round(cbind(vl0$lims[, 1], limsd, vl0$lims[, 2], pct), rou)
    limits <- cbind(vl0$lims[, 1], limsd, vl0$lims[, 2], pct)
    dimnames(limits) <- list(alpha, c("bca", "jacksd", "std", "pct"))
    ##stats <- round(rbind(stats0, statsd), rou)
    stats <- rbind(stats0, statsd)
    dimnames(stats) <- list(c("est", "jsd"), c("theta", "sdboot", "z0", "a", "sdjack"))
    vl <- list(call = call, lims = limits, stats = stats, B.mean = B.mean, seed = seed)
    if (ttind == 0) {
        ##vl$ustats <- round(ustats, rou)
        vl$ustats <- ustats
    }
    ## if (sw == 5) {
    ##     vl$tt <- tt
    ##     return(vl)
    ## }
    bcaboot.return(vl)
}
