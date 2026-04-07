## Version of June 16, 2018
##
## Nonparametric BCa via regression on bootstrap count vectors.
## Estimates the acceleration 'a' by regressing bootstrap statistics on
## count vectors near uniform (1,...,1), as opposed to bcajack which uses
## classical jackknife deletion. Allows pre-computed Blist for distributed
## computation.

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
bcajack2 <- function(x, B, func, ..., m = NULL, mr, pct = 0.333, K = 2, J = 12,
                     alpha = c(0.025, 0.05, 0.1, 0.16),
                     verbose = TRUE) {

    call <- match.call()
    if (!missing(mr))
        warning("'mr' argument is currently unused in bcajack2")

    ## Save rng state
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        stats::runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    ## Inner function: compute BCa quantities from count matrix Y and
    ## bootstrap replicates tt. Estimates influence via local regression
    ## on count vectors nearest to uniform (1,...,1). EN20 appendix.
    qbca2 <- function(Y, tt, t0, alpha, pct) {
        m <- ncol(Y)
        B <- nrow(Y)

        ra <- regression_accel(Y, tt, t0, pct)
        reg_infl <- ra$reg_infl
        a <- ra$a
        sdjack <- ra$sdjack

        s <- mean(tt)
        B.mean <- c(B, s)

        zalpha <- stats::qnorm(alpha)
        n_alpha <- length(alpha)
        ## cov_infl = bootstrap covariance-based influence (m * Cov(tt, Y))
        cov_infl <- m * as.vector(stats::cov(tt, Y))
        ## Variance-stabilized influence: 2*regression - bootstrap_cov. EN20 Sec 3
        adj_infl <- 2 * reg_infl - cov_infl
        ustats <- compute_ustats(t0, tt, adj_infl, m)

        sdboot <- stats::sd(tt)
        ## Bias-correction z0. Efron (1987) Sec 2
        z0 <- stats::qnorm(sum(tt < t0)/B)

        ## BCa percentile formula. Efron (1987) Sec 2
        lims <- compute_bca_limits(z0, a, zalpha, tt)$limits
        ## Standard (normal-theory) limits for comparison
        standard <- t0 + sdboot * stats::qnorm(alpha)
        lims <- cbind(lims, NA_real_, standard, NA_real_)
        dimnames(lims) <- list(alpha, c("bca", "jacksd", "std", "pct"))
        stats <- c(t0, sdboot, z0, a, sdjack)
        names(stats) <- c("theta", "sdboot", "z0", "a", "sdjack")
        list(lims = lims, stats = stats, B.mean = B.mean, ustats = ustats)
    }

    alpha <- expand_alpha(alpha)

    boot_data <- if (is.list(B)) B else NULL
    bs <- bootstrap_resample(x, B, func, ..., m = m, verbose = verbose,
                             boot_data = boot_data)
    Y <- bs$Y; tt <- bs$tt; t0 <- bs$t0; B <- bs$B
    result0 <- qbca2(Y, tt, t0, alpha = alpha, pct = pct)
    result0$call <- call
    result0$seed <- seed

    if (K == 0)
        bcaboot.return(result0)

    n_alpha <- length(alpha)
    Pct <- numeric(n_alpha)
    for (i in seq_len(n_alpha)) Pct[i] <- sum(tt <= result0$lims[i, 1])/B
    Stand <- result0$stats[1] + result0$stats[2] * stats::qnorm(alpha)
    Limsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)

    ## Internal SE via delete-d jackknife of the bootstrap.
    ## Split B replications into J groups, recompute qbca2 leaving one group
    ## out at a time. Repeat K times, average. Scale: (J-1)/sqrt(J).
    Limbcsd <- matrix(0, length(alpha), K)
    Statsd <- matrix(0, 5, K)
    for (k in seq_len(K)) {
        remainder <- B%%J
        fold_idx <- sample(x = B, size = B - remainder)
        fold_idx <- matrix(fold_idx, ncol = J)
        limbc <- limst <- matrix(0, length(alpha), J)
        stats <- matrix(0, 5, J)

        for (j in seq_len(J)) {
            iij <- c(fold_idx[, -j])
            Yj <- Y[iij, ]
            ttj <- tt[iij]
            fold_result <- qbca2(Yj, ttj, t0, alpha, pct)
            limbc[, j] <- fold_result$lims[, 1]
            limst[, j] <- fold_result$lims[, "std"]
            stats[, j] <- fold_result$stats
        }

        ## Delete-d jackknife scale correction: (J-1)/sqrt(J)
        Limbcsd[, k] <- apply(limbc, 1, sd) * (J - 1)/sqrt(J)
        Statsd[, k] <- apply(stats, 1, sd) * (J - 1)/sqrt(J)
    }
    limsd <- rowMeans(Limbcsd, 1)
    statsd <- rowMeans(Statsd, 1)
    limits <- cbind(result0$lims[, "bca"], limsd, result0$lims[, "std"], Pct)
    dimnames(limits) <- list(alpha, c("bca", "jacksd", "std", "pct"))
    stats <- rbind(result0$stats, statsd)
    ustats <- result0$ustats
    B.mean <- c(B, mean(tt))
    dimnames(stats) <- list(c("est", "jsd"), c("theta", "sdboot", "z0", "a", "sdjack"))
    result <- list(call = call, lims = limits, stats = stats, B.mean = B.mean, ustats = ustats,
                seed = seed)
    bcaboot.return(result)
}
