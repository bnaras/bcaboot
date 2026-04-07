#' Compute nonparametric BCa bootstrap confidence intervals
#'
#' This function computes nonparametric BCa confidence intervals using
#' regression on bootstrap count vectors (same approach as [bcajack2])
#' plus a generalized BCa (gbca) diagnostic that tests adequacy of the
#' BCa approximation via warped normal comparison. Efron (1987) Sec 7.
#'
#' @param B number of bootstrap replications, or a list with components
#'   `Y` (count matrix), `tt` (bootstrap replicates), `t0` (original estimate)
#' @param x data matrix of dimension \eqn{n\times p}, rows assumed mutually independent
#' @param func R function, where \eqn{func(x, \ldots)} is the statistic of interest
#' @param ... additional arguments for `func`
#' @param m number of units; if \eqn{m<n} then original units collected into
#'   groups each of size \eqn{n/m}; useful if \eqn{n} is very large
#' @param pct proportion of nearby count vectors used in finding gradient
#' @param K number of jackknife repetitions for internal standard error estimation
#' @param J number of jackknife folds per repetition
#' @param alpha percentiles of coverage probabilities for bootstrap intervals
#' @param verbose logical for verbose progress messages
#' @return a named list of class `bcaboot` containing:
#'
#' * __call__: the matched call
#'
#' * __lims__ : BCa confidence limits (`bca`), jackknife internal
#'     standard errors (`jacksd`), standard limits (`std`), and
#'     bootstrap percentiles (`pct`)
#'
#' * __stats__ : Estimates and jackknife Monte Carlo errors for
#'     `theta`, `sdboot`, `sdjack`, `z0`, `a`, plus big-A acceleration
#'     and SE stability (`A`, `sese`)
#'
#' * __B.mean__ : bootstrap sample size B and mean of replications
#'
#' * __equiv__ : gbca coverage estimates showing `alpha` vs `alphatilda`
#'     (coverage-adjusted alpha); if BCa is exact these are equal.
#'     Also jackknife Monte Carlo standard deviations.
#'
#' * __seed__ : the random number state for reproducibility
#'
#' @importFrom stats approx smooth.spline
#' @export
bcanon <- function (B, x, func, ..., m = nrow(x), pct = .333, K = 2, J = 12,
                    alpha = c(0.025, 0.05, 0.1, 0.16),
                    verbose = TRUE) {

    call <- match.call()
    alpha <- expand_alpha(alpha)

    ## Save rng state
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        stats::runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    ## Inner function: compute BCa quantities plus generalized BCa (gbca)
    ## diagnostic. Uses regression on nearby count vectors (same approach as
    ## bcajack2's qbca2) plus the gbca warping from Efron (1987) Sec 7.
    ## Accesses alpha, pct from enclosing scope.
    qbca2 <- function(Y, tt, t0) {
        m <- ncol(Y)
        B <- nrow(Y)

        ra <- regression_accel(Y, tt, t0, pct)
        reg_infl <- ra$reg_infl
        a <- ra$a
        sdjack <- ra$sdjack
        local_dir <- ra$local_dir

        s <- mean(tt)
        B.mean <- c(B, s)
        zalpha <- qnorm(alpha)
        n_alpha <- length(alpha)
        sdboot <- sd(tt)
        ## Bias-correction z0. Efron (1987) Sec 2
        z0 <- qnorm(sum(tt < t0)/B)

        ## Big-A raw acceleration via delta method. Efron (1987) Sec 7.
        std_tt <- (tt - mean(tt))/sdboot
        A <- 0.5 * (sdboot/sdjack) * sum(std_tt^2 * local_dir)/B
        ## rms = RMS of covariance between centered counts and squared standardized
        ## bootstrap replicates; measures stability of SE across data perturbations
        Del <- cov(YY, std_tt^2)/2
        rms <- sqrt((m/(m - 1)) * sum(Del^2))
        sese <- c(A, rms)
        names(sese) <- c("A", "rms")

        ## BCa percentile formula. Efron (1987) Sec 2
        lims <- compute_bca_limits(z0, a, zalpha, tt)$limits
        ## Standard (normal-theory) limits for comparison
        standard <- t0 + sdboot * qnorm(alpha)
        lims <- cbind(lims, standard)
        dimnames(lims) <- list(alpha, c("bca", "std"))
        stats <- c(t0, sdboot, z0, a, sdjack)
        names(stats) <- c("thet", "sdboot", "z0", "a", "sdjack")
        result <- list(lims = lims, stats = stats, B.mean = B.mean, sese = sese,
                       t0 = t0, tt = tt, local_dir = local_dir)

        ## ---------------------------------------------------------------
        ## Generalized BCa (gbca) diagnostic. Efron (1987) Sec 7.
        ## Tests adequacy of the BCa approximation by comparing with a
        ## warped normal model. If BCa is exact, the diagnostic D(z)
        ## should be constant = 1.
        ## ---------------------------------------------------------------

        ## Grid of standard normal quantiles for diagnostic evaluation
        zz <- seq(-2.5, 2.5, 0.05)
        cdf_grid <- pnorm(zz)
        n_grid <- length(cdf_grid)  # 101 points; index 51 is the center (z=0)
        tta <- as.vector(quantile(tt, cdf_grid))
        ## ecdf_vals[i] = empirical CDF of local_dir up to the
        ## bootstrap quantile at level cdf_grid[i]
        ecdf_vals <- numeric(n_grid)
        for (i in seq_len(n_grid)) ecdf_vals[i] <- sum(local_dir[tt <= tta[i]])/B
        ## H = ratio of empirical CDF to normal density;
        ## diag_raw = H normalized at center (z=0, index 51). D(z)=1 means BCa exact.
        H <- ecdf_vals/dnorm(zz)
        diag_raw <- H/H[51]
        diag_smooth <- smooth.spline(zz, diag_raw, df = 5)$y

        Dmatrix <- cbind(zz, diag_raw, diag_smooth)
        z0 <- stats[3]
        a0 <- stats[4]
        ## Scaled acceleration for warping: scaled_a = a0/(1-a0*z0)
        scaled_a <- a0/(1 - a0 * z0)
        i0 <- 51  # center index in zz grid

        ## Warping transformation: integrates 1/diag_smooth (the inverse diagnostic)
        ## to build a nonlinear scale on which the BCa approximation would
        ## be exact. When |scaled_a| is near 0, ww ~ I (identity-like).
        I <- (cumsum(1/diag_smooth) - 0.5/diag_smooth) * diff(zz)[1]
        I <- I - I[i0]
        if (abs(scaled_a) < 1e-06)
            ww <- I else ww <- (exp(scaled_a * I) - 1)/scaled_a

        cdf_grid <- pnorm(zz)

        ## BCa transformation on original scale
        Z <- z0 + (z0 + zz)/(1 - a0 * (z0 + zz))
        bet <- pnorm(Z)

        ## Warping functions: wfun maps z -> warped scale,
        ## winv is the inverse, Phitil = Phi(winv(w)) = CDF on warped scale
        wfun <- function(z) {
            approx(zz, ww, z, rule = 2, ties = mean)$y
        }
        winv <- function(w) {
            approx(ww, zz, w, rule = 2, ties = mean)$y
        }
        Phitil <- function(w) {
            pnorm(winv(w))
        }

        ## Apply BCa formula on warped scale to get gbca coverage
        zt0 <- wfun(z0)
        zzt <- -wfun(rev(zz))
        Zt <- zt0 + (zt0 + zzt)/(1 - a0 * (zt0 + zzt))
        bett <- Phitil(Zt)

        ## equiv: alpha-tilde = gbca-adjusted coverage corresponding to
        ## nominal BCa coverage alpha. If BCa is exact, alpha-tilde = alpha.
        atil <- approx(bett, cdf_grid, bet, rule = 2, ties = mean)$y
        al <- alpha
        altil <- approx(cdf_grid, atil, al, rule = 2, ties = mean)$y
        equiv <- rbind(al, altil)
        dimnames(equiv)[[2]] <- rep(" ", length(alpha))
        dimnames(equiv)[[1]] <- c("alpha", "alphatil")
        result$equiv <- equiv
        result$Dmatrix <- cbind(zz, diag_smooth)
        return(result)
    }

    boot_data <- if (is.list(B)) B else NULL
    bs <- bootstrap_resample(x, B, func, ..., m = m, verbose = verbose,
                             boot_data = boot_data)
    Y <- bs$Y; tt <- bs$tt; t0 <- bs$t0; B <- bs$B
    result0 <- qbca2(Y, tt, t0)
    result0$call <- call

    if (K == 0)
        return(result0)
    equiv <- result0$equiv
    n_alpha <- length(alpha)
    Pct <- numeric(n_alpha)
    for (i in seq_len(n_alpha)) Pct[i] <- sum(tt <= result0$lims[i, 1])/B
    Stand <- result0$stats[1] + result0$stats[2] * qnorm(alpha)
    Eqsd <- matrix(0, n_alpha, K)

    ## Internal SE via delete-d jackknife of the bootstrap.
    ## Split B replications into J groups, recompute qbca2 leaving one group
    ## out at a time. Repeat K times, average. Scale: (J-1)/sqrt(J).
    Limbcsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)
    sesed <- matrix(0, 2, K)
    for (k in seq_len(K)) {
        remainder <- B%%J
        fold_idx <- sample(B, B - remainder)
        fold_idx <- matrix(fold_idx, ncol = J)
        limbc <- limst <- matrix(0, length(alpha), J)
        stats <- matrix(0, 5, J)
        armsd <- matrix(0, 2, J)
        eqsd <- matrix(0, n_alpha, J)
        for (j in seq_len(J)) {
            iij <- c(fold_idx[, -j])
            Yj <- Y[iij, ]
            ttj <- tt[iij]
            fold_result <- qbca2(Yj, ttj, t0)
            limbc[, j] <- fold_result$lims[, 1]
            limst[, j] <- fold_result$lims[, 2]
            stats[, j] <- fold_result$stats
            armsd[, j] <- fold_result$sese
            eqsd[, j] <- fold_result$equiv[2, ]
        }
        ## Delete-d jackknife scale correction: (J-1)/sqrt(J)
        sesed[, k] <- apply(armsd, 1, sd) * (J - 1)/sqrt(J)
        Eqsd[, k] <- apply(eqsd, 1, sd) * (J - 1)/sqrt(J)
        Limbcsd[, k] <- apply(limbc, 1, sd) * (J - 1)/sqrt(J)
        Statsd[, k] <- apply(stats, 1, sd) * (J - 1)/sqrt(J)
    }
    Armjsd <- rowMeans(sesed, 1)
    Armat <- rbind(result0$sese, Armjsd)
    dimnames(Armat) <- list(c("est", "jsd"), c("A", "sese"))
    eqjsd <- rowMeans(Eqsd, 1)
    equiv <- rbind(equiv, eqjsd)
    dimnames(equiv)[[1]] <- c("alpha", "alphatilda", "jacksd")
    limsd <- rowMeans(Limbcsd, 1)
    statsd <- rowMeans(Statsd, 1)
    limits <- cbind(result0$lims[, 1], limsd, result0$lims[, 2], Pct)
    dimnames(limits) <- list(alpha, c("bca", "jacksd", "std", "pct"))
    stats <- rbind(result0$stats, statsd)
    B.mean <- c(B, mean(tt))
    dimnames(stats) <- list(c("estimate", "jacksd"), c("theta", "sdboot", "z0", "a", "sdjack"))
    stats <- cbind(stats[, c(1, 2, 5, 3, 4)], Armat)
    result <- list(call = call, lims = limits, stats = stats, B.mean = B.mean,
                   equiv = equiv, seed = seed)
    bcaboot.return(result)
}
