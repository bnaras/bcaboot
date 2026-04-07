#' Compute nonparametric bca bootstrap confidence intervals
#' @param B number of bootstrap replications
#' @param x data matrix of dimension \eqn{n\times p}, rows assumed mutually independent
#' @param func R function, where \eqn{func(x, \ldots)} is the statistic of interest
#' @param m number of units; if \eqn{m<n} then original units collected into groups each of size \eqn{n/m}; useful if \eqn{n} is very large
#' @param pct proportion of original units used in finding gradient
#' @param K controls jackknife standard error calculations (how exactly?)
#' @param J controls jackknife standard error calculations (how exactly?)
#' @param alpha percentiles of coverage probabilities for bootstrap intervals
#' @param rou number of digits rounded to
#' @param catj counts number of bootstrap replictions; a value of 0 suppresses count
#' @param sw if nonzero, adds the replicates `tt` and `dd` to output; refer to the returned value below
#' @return a named list containing:
#'
#' * __call__: the call
#'
#' * __lims__ : Bca confidence limits (first column) and the standard
#'     limits (third column). The second column, `jacksd`, are the
#'     jackknife estimates of Monte Carlo error; `pct`, the third
#'     column are the proportion of the replicates `tt` less than each
#'     `bcalim` value
#'
#' * __stats__ : Estimates and their jackknife Monte Carlo errors:
#'     `theta` = \eqn{\hat{\theta}}; `sd`, the bootstrap standard
#'     deviation for \eqn{\hat{\theta}}; `a0` the acceleration
#'     estimate; `A` the big-A measure of raw acceleration; `sese` of
#'     standard error stablity
#'
#' * __equiv__ : shows the `gbca` coverage estimates
#'     $\eqn{\tilde{\alpha}}$ corresponding to the nomimal bca
#'     coverages \eqn{\alpha}; also their jacknife montecarlo standard
#'     deviations.
#' @export
bcanon <- function (B, x, func, ..., m = nrow(x), pct = .333, K = 2, J = 12,
                    alpha = c(0.025, 0.05, 0.1, 0.16),
                    rou=3, catj=500, sw = 0) {
  # B=number of bootstrap replications x =nxp data matrix, rows assumed mutually independent
  # func = R function, func(x) the statistic of interest m=number of units; if m<n then
  # original units collected into groups each of size n/m; useful if n is very large pct=
  # proportion of original units used in finding gradient K,J controls jackknife standard
  # error calculations alpha= percentiles of coverage probabilites for bootstrap intervals rou
  # = number of digits rounded to catj counts number of bootstrap replictions; catj=0
  # suppresses count sw=1 adds tt and dd to output
  
  # output: the call (includes B, the number of bootstrap replications) lims gives bootstrap
  # and standard interval limts, and the jackknife montecarlo standard deviations for the boot
  # limits pct if the percentiles of the B boots giving the boot limits stats gives estimates
  # of theta, its bootstrap sd, z0, and a, measures A and sese of standard error stablity, as
  # well as jackknife monte carlo sds for all of these equiv shows the gbca coverage estimates
  # alphatildacorresponding to the nomimal bca coverages alpha; also their jacknife montecarlo
  # sds
  
  call <- match.call()
  alpha <- alpha[alpha < 0.5]
  alpha <- c(alpha, 0.5, rev(1 - alpha))
  
  ## Save rng state
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    stats::runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  ##qbca2 <- function(Y, tt, t0, alpha = alpha, pct = pct, rou = rou, sw = sw) {
  qbca2 <- function(Y, tt, t0) {    
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
        sdboot <- sd(tt)
        sdjack <- sqrt(sum(ty.^2))/(m - 1)
        z0 <- qnorm(sum(tt < t0)/B)
        ## big-A acceleration
        YY <- scale(Y, T, F)
        dd <- as.vector(YY %*% ty.)/(m - 1)
        dd <- dd/sd(dd)
        DD <- (tt - mean(tt))/sdboot
        A <- 0.5 * (sdboot/sdjack) * sum(DD^2 * dd)/B
        Del <- cov(YY, DD^2)/2
        rms <- sqrt((m/(m - 1)) * sum(Del^2))
        sese <- c(A, rms)
        names(sese) <- c("A", "rms")
        iles <- pnorm(z0 + (z0 + zalpha)/(1 - a * (z0 + zalpha)))
        ooo <- trunc(iles * B)
        ooo <- pmin(pmax(ooo, 1), B)
        lims <- sort(tt)[ooo]
        standard <- t0 + sdboot * qnorm(alpha)
        lims <- round(cbind(lims, standard), rou)
        dimnames(lims) <- list(alpha, c("bcalims", "std"))
        stats <- round(c(t0, sdboot, z0, a, sdjack), rou)
        names(stats) <- c("thet", "sdboot", "z0", "a", "sdjack")
        vl <- list(lims = lims, stats = stats, B.mean = B.mean, sese = sese, t0 = t0, tt = tt, dd = dd)

        zz <- seq(-2.5, 2.5, 0.05)
        aa <- pnorm(zz)
        na <- length(aa)
        tta <- as.vector(quantile(tt, aa))
        F. <- rep(0, na)
        for (i in 1:na) F.[i] <- sum(dd[tt <= tta[i]])/B
        H <- F./dnorm(zz)
        D <- H/H[51]
        Dsm <- smooth.spline(zz, D, df = 5)$y
        DD <- cbind(D, Dsm)
        dimnames(DD)[[2]] <- c("Diagnostic", "Smoothed")

        D9 <- DD[, 2]
        Dmatrix <- cbind(zz, DD)
        z0 <- stats[3]
        a0 <- stats[4]
        ep <- a0/(1 - a0 * z0)
        i0 <- 51

        I <- (cumsum(1/D9) - 0.5/D9) * diff(zz)[1]
        I <- I - I[i0]
        if (abs(ep) < 1e-06)
            ww <- I else ww <- (exp(ep * I) - 1)/ep

        z0 <- round(z0, 3)
        a0 <- round(a0, 3)
        aa <- pnorm(zz)
        aabar <- 1 - aa

        Z <- z0 + (z0 + zz)/(1 - a0 * (z0 + zz))
        bet <- pnorm(Z)

        wfun <- function(z) {
            approx(zz, ww, z, rule = 2, ties = mean)$y
        }
        winv <- function(w) {
            approx(ww, zz, w, rule = 2, ties = mean)$y
        }
        Phitil <- function(w) {
            pnorm(winv(w))
        }

        zt0 <- wfun(z0)
        zzt <- -wfun(rev(zz))
        Zt <- zt0 + (zt0 + zzt)/(1 - a0 * (zt0 + zzt))
        bett <- Phitil(Zt)

        atil <- approx(bett, aa, bet, rule = 2, ties = mean)$y
        al <- alpha
        altil <- approx(aa, atil, al, rule = 2, ties = mean)$y
        equiv <- round(rbind(al, altil), 5)
        dimnames(equiv)[[2]] <- rep(" ", length(alpha))
        dimnames(equiv)[[1]] <- c("alpha", "alphatil")
        vl$equiv <- equiv
        vl$Dmatrix <- cbind(zz, D9)
        return(vl)
    }

    if (is.list(B)) {
        Y <- B$Y
        tt <- B$tt
        t0 <- B$t0
        B <- length(tt)
        vl0 <- qbca2(Y, tt, t0)
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
            vl0 <- qbca2(Y, tt, t0)
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
            vl0 <- qbca2(Y, tt, t0)
        }
    }
    vl0$call <- call

    if (K == 0)
        return(vl0)
    equiv <- vl0$equiv
    nal <- length(alpha)
    Pct <- rep(0, nal)
    for (i in 1:nal) Pct[i] <- round(sum(tt <= vl0$lims[i, 1])/B, 3)
    Stand <- vl0$stats[1] + vl0$stats[2] * qnorm(alpha)
    Limsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)
    Eqsd <- matrix(0, nal, K)  ###**

    Limbcsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)
    sesed <- matrix(0, 2, K)
    for (k in 1:K) {
        rr <- B%%J
        II <- sample(B, B - rr)
        II <- matrix(II, ncol = J)
        limbc <- limst <- matrix(0, length(alpha), J)
        stats <- matrix(0, 5, J)
        armsd <- matrix(0, 2, J)
        eqsd <- matrix(0, nal, J)
        for (j in 1:J) {
            iij <- c(II[, -j])
            Yj <- Y[iij, ]
            ttj <- tt[iij]
            vlj <- qbca2(Yj, ttj, t0)
            limbc[, j] <- vlj$lims[, 1]
            limst[, j] <- vlj$lims[, 2]
            stats[, j] <- vlj$stats
            armsd[, j] <- vlj$sese
            eqsd[, j] <- vlj$equiv[2, ]
        }
        sesed[, k] <- apply(armsd, 1, sd) * (J - 1)/sqrt(J)
        Eqsd[, k] <- apply(eqsd, 1, sd) * (J - 1)/sqrt(J)
        Limbcsd[, k] <- apply(limbc, 1, sd) * (J - 1)/sqrt(J)
        Statsd[, k] <- apply(stats, 1, sd) * (J - 1)/sqrt(J)
        if (catj >= 0)
            cat("{", k, "}", sep = "")
    }
    Armjsd <- rowMeans(sesed, 1)
    Armat <- rbind(vl0$sese, Armjsd)
    Armat <- round(Armat, 3)
    dimnames(Armat) <- list(c("est", "jsd"), c("A", "sese"))
    eqjsd <- rowMeans(Eqsd, 1)
    equiv <- round(rbind(equiv, eqjsd), rou)
    dimnames(equiv)[[1]] <- c("alpha", "alphatilda", "jacksd")
    limsd <- rowMeans(Limbcsd, 1)
    statsd <- rowMeans(Statsd, 1)
    limits <- round(cbind(vl0$lims[, 1], limsd, vl0$lims[, 2], Pct), rou)
    dimnames(limits) <- list(alpha, c("bcalims", "jacksd", "std", "Pct"))
    stats <- round(rbind(vl0$stats, statsd), rou)
    B.mean <- c(B, round(mean(tt), rou))
    dimnames(stats) <- list(c("estimate", "jacksd"), c("thet", "sdboot", "z0", "a0", "sdbar"))
    stats <- cbind(stats[, c(1, 2, 5, 3, 4)], Armat)
    vll <- list(call = call, lims = limits, stats = stats, equiv = equiv)
    if (sw == 1) {
        vll$tt <- tt
        vll$dd <- vl0$dd
    }
    return(vll)
}
