## Version of June 1, 2018

#'
#' @title Nonparametric bias-corrected and accelerated bootstrap
#'     confidence limits
#'
#' @description \code{bcajack} computes nonparametric confidence
#'     intervals for bootstrap estimates. For reproducibility, save or
#'     set the random number state before calling this routine.
#'
#' @details
#' Bootstrap confidence intervals depend on three elements:
#' \itemize{
#' \item the cdf of the B bootstrap replications \code{t[i]*, i=1...B}
#' \item the bias-correction number \code{z0=qnorm(#{t[i]* < t0}/B )}
#' where \code{t0=func(x)} is the original estimate
#' \item the acceleration number \eqn{a} that measures the rate of
#' change in \code{sd{t0}} as x, the data changes.
#' }
#'
#' The first two of these depend only on the bootstrap distribution,
#' and not how it is generated: parametrically or
#' non-parametrically. Program bcajack can be used in a hybrid fashion
#' in which the vector \code{tt} of B bootstrap replications is first
#' generated from a parametric model.
#'
#' So, in the diabetes example below, we might first draw bootstrap
#' samples \eqn{y^* \sim N(X\hat{\beta}, \hat{\sigma}^2I)} where
#' \eqn{\hat{beta}} and \eqn{\hat{sigma}} were obtained from
#' \code{lm(y~X)}; each \eqn{y^*} would then provide a bootstrap
#' replication \code{t* = rfun(cbind(X,y*))}.  Then we could get bca
#' intervals from \code{bcajack(Xy, tt, rfun ....)} with \code{tt},
#' the vector of B \code{t*} values. The only difference from a full
#' parametric bca analysis would lie in the nonparametric estimation
#' of \code{a}, often a negligible error.
#'
#' @param x an nxp data matrix, rows are observed p-vectors, assumed
#'     to be independently sampled from target population. If `p` is 1
#'     then x can be a vector.
#' @param B number of bootstrap replications. 'B' can also be a vector
#'     of B bootstrap replications of the estimated parameter of
#'     interest, computed separately.
#' @param func function \eqn{\hat{\theta}}=func(x) computing estimate of the
#'     parameter of interest; func(x) should return a real value for
#'     any n' x p matrix x', n' not necessarily equal to n.
#' @param ... additional arguments for func.
#' @param m integer m <= n; collects the n rows of x into m groups to
#'     speed up the jackknife calculations for estimating the
#'     acceleration value 'a'; typically m=20 or 40; does not have to
#'     exactly divide n.
#' @param mr if m < n then mr repetions of the randomly grouped
#'     jackknife calculations are averaged.
#' @param K If K > 0, bcajack also returns estimates of 'internal
#'     standard error', that is, of the variability due to stopping at
#'     B bootstrap replications rather than going on to
#'     infinity. These are obtained from a second type of jackknifing,
#'     taking an average of K separate jackknife estimates, each
#'     randomly splitting the B bootstrap replications into J groups.
#' @param J the number of groups into which the bootstrap replications are split
#' @param alpha percentiles desired for the bca confidence limits.
#' @param rou rounding parameter for the output.
#' @param verbose logical for verbose messages
#' @param sw switch that controls output, eg sw=5 returns the B
#'     bootstrap replications as well as the bca output.
#' @return a named list of several items:
#'
#' \describe{
#' \item{lims}{first column shows the estimated bca confidence limits
#' at the requested alpha percentiles. These can be compared with the
#' standard limits \code{thetahat+sdboot*z[alpha]}, third column. The second
#' column "jacksd" gives the internal standard errors for the bca
#' limits, quite small in the example. Column 4, "pct", gives the
#' percentiles of the ordered B bootstrap replications corresponding
#' to the bca limits, eg the 897th largest replication equalling the
#' .975 bca limit .557.}
#' \item{stats}{ top line of stats shows 5 estimates: thet is func(x),
#' original point estimate of the parameter of interest; sdboot is its
#' bootstrap estimate of standard error; z0 is the bca bias correction
#' value, in this case quite negative; a is the "acceleration", a
#' component of the bca limits (nearly zero here); sdjack is the
#' jackknife estimate of standard error for thet. Bottom line gives
#' the internal standard errors for the five quantities above. This is
#' substantial for z0 above.}
#' \item{B.mean}{bootstrap sample size B, and the mean of the B
#' bootstrap replications theathat*.}
#' }
#'
#' @references DiCiccio T and Efron B (1996). Bootstrap confidence
#'     intervals. Statistical Science 11, 189-228
#' @references Efron B (1987). Better bootstrap confidence
#'     intervals. JASA 82 171-200
#'
#'
#' @import lars
#' @export
#' @examples
#' data(diabetes, package = "lars")
#' Xy <- cbind(diabetes$x, diabetes$y)
#' rfun <- function(Xy) {
#'   y <- Xy[, 11]
#'   X <- Xy[, 1:10]
#'   summary(lm(y~X) )$adj.r.squared
#' }
#' bcajack(x = Xy, B = 1000, func = rfun, m = 40)
#'
bcajack <- function(x, B, func, ..., m = nrow(x), mr = 5, K = 2, J = 10, alpha = c(0.025,
    0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95, 0.975),  verbose = TRUE, sw = 0) {
    ## x is nxp data matrix, func is statistic thetahat=func(x) can enter #bootsize B
    ## for bootsim vector tt (which is calculated)

    call <- match.call()
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
    u <- rep(0, m)
    m1 <- sqrt(m * (m - 1))

    if (m == n) {
        for (i in 1:n) {
            u[i] <- func(x[-i, ], ...)
        }
        t. <- (mean(u) - u) * (m - 1)
        a <- (1/6) * sum(t.^3)/(sum(t.^2))^1.5
        sdjack <- sqrt(sum(t.^2))/m1
    }

    if (m < n) {
        aa <- ssj <- rep(0, mr)
        for (k in 1:mr) {
            r <- n%%m
            Imat <- matrix(sample(1:n, n - r), m)
            Iout <- setdiff(1:n, Imat)
            if (sw == 3)
                return(list(Imat = Imat, Iout = Iout))
            for (j in 1:m) {
                Ij <- setdiff(1:m, j)
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
    if (sw == 2)
        return(list(t. = t., sdjack = sdjack))

    if (ttind == 0) {
        tY. <- Y. <- rep(0, n)
        if (verbose) pb <- txtProgressBar(min = 0, max = B, style = 3)
        for (j in seq_len(B)) {
            ij <- sample(n, n, T)
            Yj <- table(c(ij, 1:n)) - 1
            tt[j] <- func(x[ij, ], ...)
            tY. <- tY. + tt[j] * Yj
            Y. <- Y. + Yj
            if (verbose) setTxtProgressBar(pb, j)
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
    zalpha <- qnorm(alpha)
    nal <- length(alpha)

    sdboot0 <- sd(tt)  # sdd=sd(dd)
    z00 <- qnorm(sum(tt < t0)/B)

    iles <- pnorm(z00 + (z00 + zalpha)/(1 - a * (z00 + zalpha)))
    ooo <- trunc(iles * B)
    ooo <- pmin(pmax(ooo, 1), B)
    lims0 <- sort(tt)[ooo]
    standard <- t0 + sdboot0 * qnorm(alpha)
    ## lims0 <- round(cbind(lims0, standard), rou)
    lims0 <- cbind(lims0, standard)
    dimnames(lims0) <- list(alpha, c("bcalims", "standard"))
    ## stats0 <- round(c(t0, sdboot0, z00, a, sdjack), rou)
    stats0 <- c(t0, sdboot0, z00, a, sdjack)
    names(stats0) <- c("thet", "sdboot", "z0", "a", "sdjack")
    vl0 <- list(lims0 = lims0, stats0 = stats0, B.mean = B.mean, call = call)
    if (K == 0)
        return(vl0)

    pct <- rep(0, nal)
    ##for (i in 1:nal) pct[i] <- round(sum(tt <= lims0[i, 1])/B, 3)
    for (i in 1:nal) pct[i] <- sum(tt <= lims0[i, 1])/B
    Stand <- vl0$stats[1] + vl0$stats[2] * qnorm(alpha)
    Limsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)

    for (k in 1:K) {
        II <- sample(B, B)
        II <- matrix(II, ncol = J)
        lims <- matrix(0, length(alpha), J)
        stats <- matrix(0, 5, J)
        for (j in 1:J) {
            iij <- c(II[, -j])
            ttj <- tt[iij]
            Bj <- length(ttj)
            sdboot <- sd(ttj)
            z0 <- qnorm(sum(ttj < t0)/Bj)

            iles <- pnorm(z0 + (z0 + zalpha)/(1 - a * (z0 + zalpha)))
            oo <- trunc(iles * Bj)
            oo <- pmin(pmax(oo, 1), Bj)
            li <- sort(ttj)[oo]
            standard <- t0 + sdboot * qnorm(alpha)
            ##sta <- round(c(t0, sdboot, z0, a, sdjack), rou)
            sta <- c(t0, sdboot, z0, a, sdjack)
            names(sta) <- c("thet", "sdboot", "z0", "a", "sdjack")
            lims[, j] <- li
            stats[, j] <- sta
        }
        if (sw == 4)
            return(list(lims = lims, stats = stats))
        Limsd[, k] <- apply(lims, 1, sd) * (J - 1)/sqrt(J)
        Statsd[, k] <- apply(stats, 1, sd) * (J - 1)/sqrt(J)
        ##if (verbose) cat("{", k, "}", sep = "")
    }
    limsd <- rowMeans(Limsd, 1)
    statsd <- rowMeans(Statsd, 1)
    ##limits <- round(cbind(vl0$lims[, 1], limsd, vl0$lims[, 2], pct), rou)
    limits <- cbind(vl0$lims[, 1], limsd, vl0$lims[, 2], pct)
    dimnames(limits) <- list(alpha, c("bcalims", "jacksd", "standard", "pct"))
    ##stats <- round(rbind(stats0, statsd), rou)
    stats <- rbind(stats0, statsd)
    dimnames(stats) <- list(c("est", "jsd"), c("thet", "sdboot", "z0", "a", "sdjack"))
    vl <- list(call = call, lims = limits, stats = stats, B.mean = B.mean)
    if (ttind == 0) {
        ##vl$ustats <- round(ustats, rou)
        vl$ustats <- ustats
    }
    if (sw == 5) {
        vl$tt <- tt
        return(vl)
    }
    return(vl)
}
