## Version of June 16, 2018

#' @title Nonparametric bias-corrected and accelerated bootstrap
#'     confidence limits allowing recomputations of original statistics
#'
#' @description \code{bcajack2} is a version of \code{bcajack} that allows
#'     all the recomputations of the original statistic function
#'     \code{func} to be carried out separately. This is an advantage
#'     if \code{func} is time-consuming, in which case the B
#'     replications for the nonparametric bca calculations might need
#'     to be done on a distributed basis.
#'
#' @details
#'
#' To use bcajack2 in this mode, we first compute a list \code{Blist}
#' via \code{Blist <- list(Y = Y,tt = tt,t0 = t0)}.  Here \code{tt} is
#' a vector of length B having i-th entry \code{tt[i] <- func(x[Ii,],
#' ...)}, where x is the n x p data matrix and \code{Ii <-
#' sample(n,n,T)}, a bootstrap vector of indices. Y is a B x n "count"
#' matrix, whose i-th row is the counts corresponding to Ii. For
#' example if n = 5 and \code{Ii=(2,5,2,1,4)}, then \code{Yi =
#' (1,2,0,1,1)}. Having computed \code{Blist}, \code{bcajack2} is
#' invoked as \code{bcajack2(Blist)} without need to enter the
#' function func.
#'
#' @export
bcajack2 <- function(B, x, func, ..., m = n, pct = 0.333, K = 2, J = 12,
                     alpha = c(0.025, 0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95, 0.975),
                     rou = 3, catj = 100, sw = 0) {

    call <- match.call()

    qbca2 <- function(Y, tt, t0, alpha = c(0.025, 0.05, 0.1, 0.16, 0.5, 0.84, 0.9,
        0.95, 0.975), pct = 0.333, rou = 3, sw = 0) {
        m <- ncol(Y)
        B <- nrow(Y)
        o1 <- rep(1, m)
        D <- rep(0, B)

        for (i in 1:B) {
            Yi <- Y[i, ]
            d <- 2 * Yi * log(Yi) - 2 * (Yi - 1)
            d[Yi == 0] <- 2
            D[i] <- sum(d)
        }
        Qd <- quantile(D, pct)
        ip <- (1:B)[D <= Qd]
        ty. <- as.vector(m * lm(tt[ip] ~ Y[ip, ] - 1)$coef)
        ty. <- ty. - mean(ty.)
        a <- (1/6) * sum(ty.^3)/sum(ty.^2)^1.5
        if (sw == 3)
            return(ty.)
        s <- mean(tt)
        B.mean <- c(B, s)
        zalpha <- qnorm(alpha)
        nal <- length(alpha)
        ustat <- 2 * t0 - s
        s. <- m * .v(cov(tt, Y))
        u. <- 2 * ty. - s.
        sdu <- sum(u.^2)^0.5/m
        ustats <- c(ustat, sdu)
        names(ustats) <- c("ustat", "sdu")

        sdboot <- sd(tt)
        sdjack <- sqrt(sum(ty.^2))/(m - 1)
        z0 <- qnorm(sum(tt < t0)/B)

        iles <- pnorm(z0 + (z0 + zalpha)/(1 - a * (z0 + zalpha)))
        ooo <- trunc(iles * B)
        ooo <- pmin(pmax(ooo, 1), B)
        lims <- sort(tt)[ooo]
        standard <- t0 + sdboot * qnorm(alpha)
        lims <- round(cbind(lims, standard), rou)
        dimnames(lims) <- list(alpha, c("bcalims", "standard"))
        stats <- round(c(t0, sdboot, z0, a, sdjack), rou)
        names(stats) <- c("thet", "sdboot", "z0", "a", "sdjack")
        vl <- list(lims = lims, stats = stats, B.mean = B.mean, ustats = ustats)
        return(vl)
    }

    if (is.list(B)) {
        Y <- B$Y
        tt <- B$tt
        t0 <- B$t0
        B <- length(tt)
        vl0 <- qbca2(Y, tt, t0, alpha = alpha, pct = pct, rou = rou)
    } else {
        if (is.vector(x))
            x <- as.matrix(x)
        n <- nrow(x)
        tt <- rep(0, B)
        t0 <- func(x, ...)

        if (m == n) {
            ii <- sample(1:n, n * B, T)
            ii <- matrix(ii, B)
            Y <- matrix(0, B, n)
            for (k in 1:B) {
                ik <- ii[k, ]
                tt[k] <- func(x[ik, ], ...)
                Y[k, ] <- table(c(ik, 1:n)) - 1
                if (catj > 0)
                  if (k/catj == floor(k/catj))
                    cat("{", k, "}", sep = "")
            }
            vl0 <- qbca2(Y, tt, t0, alpha = alpha, pct = pct, rou = rou)
        }

        if (m < n) {
            r <- n%%m
            Imat <- matrix(sample(1:n, n - r), m)
            Iout <- setdiff(1:n, Imat)
            ii <- sample(1:m, m * B, T)
            ii <- matrix(ii, B)
            Y <- matrix(0, B, m)
            for (k in 1:B) {
                ik <- ii[k, ]
                Ik <- c(t(Imat[ik, ]))
                Ik <- c(Ik, Iout)
                tt[k] <- func(x[Ik, ], ...)
                Y[k, ] <- table(c(ik, 1:m)) - 1
                if (catj > 0)
                  if (k/catj == floor(k/catj))
                    cat("{", k, "}", sep = "")
            }
            vl0 <- qbca2(Y, tt, t0, alpha = alpha, pct = pct, rou = rou)
        }
    }
    vl0$call <- call

    if (K == 0)
        return(vl0)

    nal <- length(alpha)
    Pct <- rep(0, nal)
    for (i in 1:nal) Pct[i] <- round(sum(tt <= vl0$lims[i, 1])/B, 3)
    Stand <- vl0$stats[1] + vl0$stats[2] * qnorm(alpha)
    Limsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)

    Limbcsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)
    for (k in 1:K) {
        rr <- B%%J
        II <- sample(B, B - rr)
        II <- matrix(II, ncol = J)
        limbc <- limst <- matrix(0, length(alpha), J)
        stats <- matrix(0, 5, J)

        for (j in 1:J) {
            iij <- c(II[, -j])
            Yj <- Y[iij, ]
            ttj <- tt[iij]
            vlj <- qbca2(Yj, ttj, t0, alpha, pct, rou)
            limbc[, j] <- vlj$lims[, 1]
            limst[, j] <- vlj$lims[, 2]
            stats[, j] <- vlj$stats
        }

        if (sw == 4)
            return(list(limbc = limbc, limst = limst, stats = stats))
        Limbcsd[, k] <- apply(limbc, 1, sd) * (J - 1)/sqrt(J)
        Statsd[, k] <- apply(stats, 1, sd) * (J - 1)/sqrt(J)
        if (catj >= 0)
            cat("{", k, "}", sep = "")
        if (sw == 6)
            return(list(Limbcsd = Limbcsd, Statsd = Statsd))
    }
    limsd <- rowMeans(Limbcsd, 1)
    statsd <- rowMeans(Statsd, 1)
    limits <- round(cbind(vl0$lims[, 1], limsd, vl0$lims[, 2], Pct), rou)
    dimnames(limits) <- list(alpha, c("bcalims", "jacksd", "standard", "Pct"))
    stats <- round(rbind(vl0$stats, statsd), rou)
    ustats <- round(vl0$ustats, rou)
    B.mean <- c(B, round(mean(tt), rou))
    dimnames(stats) <- list(c("est", "jsd"), c("thet", "sdboot", "z0", "a", "sdjack"))
    vll <- list(call = call, lims = limits, stats = stats, B.mean = B.mean, ustats = ustats)
    if (sw == 5)
        vll$tt <- tt
    return(vll)
}


