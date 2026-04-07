#' Nonparametric BCa bootstrap confidence intervals
#'
#' Computes nonparametric bias-corrected and accelerated (BCa) bootstrap
#' confidence intervals. Consolidates the functionality of the deprecated
#' `bcajack`, `bcajack2`, and `bcanon` functions.
#'
#' @param x An \eqn{n \times p} data matrix. Rows are observed
#'   \eqn{p}-vectors, assumed independently sampled. If \eqn{p = 1}
#'   then `x` can be a vector.
#' @param B Number of bootstrap replications (integer).
#' @param func Function computing the parameter of interest.
#'   `func(x, ...)` should return a scalar for any data matrix.
#' @param ... Additional arguments passed to `func`.
#' @param accel Method for estimating acceleration:
#'   `"regression"` (default) uses local regression on bootstrap count
#'   vectors nearest to uniform; `"jackknife"` uses classical
#'   delete-one (or delete-group) jackknife influence values.
#' @param conf.level Confidence levels for the intervals (e.g.,
#'   `c(0.95, 0.90, 0.80, 0.68)`). Only values in `(0, 1)` are used.
#' @param n_groups Number of groups for jackknife. Default `nrow(x)`
#'   gives delete-one jackknife. Set smaller (e.g., 20-40) to speed up
#'   computation on large datasets.
#' @param group_reps When `n_groups < nrow(x)`, number of random
#'   grouping repetitions to average over. Only used with
#'   `accel = "jackknife"`.
#' @param kl_fraction Fraction of bootstrap count vectors nearest to
#'   uniform used in the local regression. Only used with
#'   `accel = "regression"`. Default 1/3.
#' @param n_jack Number of delete-d jackknife repetitions for estimating
#'   internal (Monte Carlo) standard error. Set to 0 to skip.
#' @param jack_groups Number of groups per jackknife fold.
#' @param boot_data Optional pre-computed bootstrap data: a list with
#'   components `Y` (B x n count matrix), `tt` (bootstrap replicates),
#'   `t0` (original estimate). When provided, `B` is ignored and no
#'   resampling is performed.
#' @param verbose Logical; show progress bar during bootstrap.
#'
#' @return An object of class `"bcaboot"` with components:
#'   \describe{
#'     \item{limits}{9-row matrix of confidence limits (columns: `bca`,
#'       `jacksd`, `std`, `pct`)}
#'     \item{stats}{2-row matrix of estimates and jackknife SEs (columns:
#'       `theta`, `sdboot`, `z0`, `a`, `sdjack`)}
#'     \item{B_mean}{Bootstrap sample size and mean of replicates}
#'     \item{ustats}{Bias-corrected estimator and its SE}
#'     \item{diagnostic}{gbca diagnostic (when `accel = "regression"`)}
#'   }
#'
#' @references Efron B (1987). Better bootstrap confidence intervals.
#'   JASA 82, 171-200.
#' @references Efron B and Narasimhan B (2020). Automatic construction
#'   of bootstrap confidence intervals. JASA 115, 1895-1905.
#'
#' @examples
#' data(diabetes, package = "bcaboot")
#' Xy <- cbind(diabetes$x, diabetes$y)
#' rfun <- function(Xy) {
#'   y <- Xy[, 11]
#'   X <- Xy[, 1:10]
#'   summary(lm(y ~ X))$adj.r.squared
#' }
#' set.seed(1234)
#' result <- bca_nonpar(x = Xy, B = 1000, func = rfun, n_groups = 34,
#'                      verbose = FALSE)
#' print(result)
#'
#' @export
bca_nonpar <- function(x, B, func, ...,
                       accel = c("regression", "jackknife"),
                       conf.level = c(0.95, 0.90, 0.80, 0.68),
                       n_groups = nrow(x),
                       group_reps = 5,
                       kl_fraction = 0.333,
                       n_jack = 2,
                       jack_groups = 12,
                       boot_data = NULL,
                       verbose = TRUE) {

    call <- match.call()
    accel <- match.arg(accel)

    ## Save rng state
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        stats::runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    if (is.vector(x))
        x <- as.matrix(x)
    n <- nrow(x)

    ## Convert conf.level to alpha and expand
    alpha <- expand_alpha((1 - conf.level) / 2)
    zalpha <- stats::qnorm(alpha)
    n_alpha <- length(alpha)

    ## --- Acceleration estimation ---
    if (accel == "jackknife") {
        ## Classical jackknife acceleration
        t0 <- func(x, ...)
        jk <- jackknife_accel(x, func, ..., m = n_groups, mr = group_reps)
        a <- jk$a
        sdjack <- jk$sdjack
        primary_infl <- jk$jack_infl

        ## Bootstrap resampling (bcajack-style: accumulate weighted sums, no Y matrix)
        ttind <- !is.null(boot_data)
        if (ttind) {
            tt <- boot_data$tt
            B <- length(tt)
        } else {
            tt <- numeric(B)
            boot_wt_sum <- boot_ct_sum <- numeric(n)
            if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
            for (j in seq_len(B)) {
                ij <- sample(x = n, size = n, replace = TRUE)
                Yj <- table(c(ij, 1:n)) - 1
                tt[j] <- func(x[ij, ], ...)
                boot_wt_sum <- boot_wt_sum + tt[j] * Yj
                boot_ct_sum <- boot_ct_sum + Yj
                if (verbose) utils::setTxtProgressBar(pb, j)
            }
            if (verbose) close(pb)
            tt_mean <- mean(tt)
            boot_wt_sum <- boot_wt_sum / B
            boot_ct_sum <- boot_ct_sum / B
            cov_infl <- n * (boot_wt_sum - tt_mean * boot_ct_sum)
            adj_infl <- 2 * primary_infl - cov_infl
            ustats <- compute_ustats(t0, tt, adj_infl, n)
        }

        ## No Y matrix or local_dir for jackknife method
        Y <- NULL
        diagnostic <- NULL

    } else {
        ## Regression-based acceleration
        bs <- bootstrap_resample(x, B, func, ..., m = n_groups, verbose = verbose,
                                 boot_data = boot_data)
        Y <- bs$Y; tt <- bs$tt; t0 <- bs$t0; B <- bs$B

        ra <- regression_accel(Y, tt, t0, kl_fraction)
        a <- ra$a
        sdjack <- ra$sdjack
        primary_infl <- ra$reg_infl
        local_dir <- ra$local_dir

        ## Ustats via variance-stabilized influence
        m <- ncol(Y)
        cov_infl <- m * as.vector(stats::cov(tt, Y))
        adj_infl <- 2 * primary_infl - cov_infl
        ustats <- compute_ustats(t0, tt, adj_infl, m)

        ## gbca diagnostic (always attempted for regression method)
        diagnostic <- tryCatch(
            gbca_diagnostic(Y, tt, t0, a, sdjack, local_dir, alpha),
            error = function(e) {
                warning("gbca diagnostic failed: ", conditionMessage(e),
                        call. = FALSE)
                NULL
            }
        )
    }

    ## --- BCa limits ---
    sdboot <- stats::sd(tt)
    z0 <- stats::qnorm(sum(tt < t0) / B)
    bca_result <- compute_bca_limits(z0, a, zalpha, tt)
    bca_lims <- bca_result$limits
    standard <- t0 + sdboot * stats::qnorm(alpha)
    B_mean <- c(B, mean(tt))

    ## Assemble K=0 result
    lims0 <- cbind(bca_lims, NA_real_, standard, NA_real_)
    dimnames(lims0) <- list(alpha, c("bca", "jacksd", "std", "pct"))
    stats0 <- c(t0, sdboot, z0, a, sdjack)
    names(stats0) <- c("theta", "sdboot", "z0", "a", "sdjack")

    if (!exists("ustats", inherits = FALSE))
        ustats <- c(ustat = NA_real_, sdu = NA_real_)

    if (n_jack == 0) {
        stats_mat <- rbind(stats0, jsd = rep(NA_real_, 5))
        dimnames(stats_mat) <- list(c("est", "jsd"),
                                    c("theta", "sdboot", "z0", "a", "sdjack"))
        return(new_bcaboot(
            limits = lims0, stats = stats_mat, B_mean = B_mean,
            ustats = ustats, call = call, seed = seed,
            method = "nonpar", accel = accel,
            diagnostic = diagnostic
        ))
    }

    ## --- Internal SE via delete-d jackknife ---
    pct <- numeric(n_alpha)
    for (i in seq_len(n_alpha)) pct[i] <- sum(tt <= lims0[i, 1]) / B

    Limsd <- matrix(0, n_alpha, n_jack)
    Statsd <- matrix(0, 5, n_jack)
    J <- jack_groups

    if (accel == "jackknife") {
        ## bcajack-style: recompute BCa limits on subsets of tt
        for (k in seq_len(n_jack)) {
            fold_idx <- matrix(sample(x = B, size = B), ncol = J)
            fold_lims <- matrix(0, n_alpha, J)
            fold_stats <- matrix(0, 5, J)
            for (j in seq_len(J)) {
                iij <- c(fold_idx[, -j])
                ttj <- tt[iij]
                Bj <- length(ttj)
                sdj <- stats::sd(ttj)
                z0j <- stats::qnorm(sum(ttj < t0) / Bj)
                fold_lims[, j] <- compute_bca_limits(z0j, a, zalpha, ttj)$limits
                fold_stats[, j] <- c(t0, sdj, z0j, a, sdjack)
            }
            Limsd[, k] <- apply(fold_lims, 1, stats::sd) * (J - 1) / sqrt(J)
            Statsd[, k] <- apply(fold_stats, 1, stats::sd) * (J - 1) / sqrt(J)
        }
    } else {
        ## bcajack2/bcanon-style: recompute full qbca2 on subsets of (Y, tt)
        ## Also collect diagnostic SEs if diagnostic was computed
        has_diag <- !is.null(diagnostic)
        if (has_diag) {
            Eqsd <- matrix(0, n_alpha, n_jack)
            sesed <- matrix(0, 2, n_jack)
        }
        for (k in seq_len(n_jack)) {
            remainder <- B %% J
            fold_idx <- matrix(sample(x = B, size = B - remainder), ncol = J)
            fold_lims <- matrix(0, n_alpha, J)
            fold_stats <- matrix(0, 5, J)
            if (has_diag) {
                fold_sese <- matrix(0, 2, J)
                fold_equiv <- matrix(0, n_alpha, J)
            }
            for (j in seq_len(J)) {
                iij <- c(fold_idx[, -j])
                Yj <- Y[iij, ]
                ttj <- tt[iij]
                Bj <- length(ttj)

                raj <- regression_accel(Yj, ttj, t0, kl_fraction)
                sdj <- stats::sd(ttj)
                z0j <- stats::qnorm(sum(ttj < t0) / Bj)
                fold_lims[, j] <- compute_bca_limits(z0j, raj$a,
                                                     stats::qnorm(alpha), ttj)$limits
                fold_stats[, j] <- c(t0, sdj, z0j, raj$a, raj$sdjack)

                if (has_diag) {
                    ## Big-A and sese for this fold
                    std_ttj <- (ttj - mean(ttj)) / sdj
                    Aj <- 0.5 * (sdj / raj$sdjack) *
                        sum(std_ttj^2 * raj$local_dir) / Bj
                    YYj <- scale(Yj, center = TRUE, scale = FALSE)
                    Delj <- stats::cov(YYj, std_ttj^2) / 2
                    mj <- ncol(Yj)
                    rmsj <- sqrt((mj / (mj - 1)) * sum(Delj^2))
                    fold_sese[, j] <- c(Aj, rmsj)
                    ## gbca equiv for this fold
                    diag_j <- tryCatch(
                        gbca_diagnostic(Yj, ttj, t0, raj$a, raj$sdjack,
                                        raj$local_dir, alpha),
                        error = function(e) NULL
                    )
                    if (!is.null(diag_j))
                        fold_equiv[, j] <- diag_j$equiv[2, ]
                }
            }
            Limsd[, k] <- apply(fold_lims, 1, stats::sd) * (J - 1) / sqrt(J)
            Statsd[, k] <- apply(fold_stats, 1, stats::sd) * (J - 1) / sqrt(J)
            if (has_diag) {
                sesed[, k] <- apply(fold_sese, 1, stats::sd) * (J - 1) / sqrt(J)
                Eqsd[, k] <- apply(fold_equiv, 1, stats::sd) * (J - 1) / sqrt(J)
            }
        }
        ## Augment diagnostic with jackknife SEs
        if (has_diag) {
            Armjsd <- rowMeans(sesed)
            sese_mat <- rbind(diagnostic$sese, Armjsd)
            dimnames(sese_mat) <- list(c("est", "jsd"), c("A", "sese"))
            diagnostic$sese <- sese_mat

            eqjsd <- rowMeans(Eqsd)
            diagnostic$equiv <- rbind(diagnostic$equiv, jacksd = eqjsd)
        }
    }

    limsd <- rowMeans(Limsd)
    statsd <- rowMeans(Statsd)

    ## Assemble final limits and stats with jackknife SEs
    limits <- cbind(lims0[, "bca"], limsd, lims0[, "std"], pct)
    dimnames(limits) <- list(alpha, c("bca", "jacksd", "std", "pct"))

    stats_mat <- rbind(stats0, statsd)
    dimnames(stats_mat) <- list(c("est", "jsd"),
                                c("theta", "sdboot", "z0", "a", "sdjack"))

    new_bcaboot(
        limits = limits, stats = stats_mat, B_mean = B_mean,
        ustats = ustats, call = call, seed = seed,
        method = "nonpar", accel = accel,
        diagnostic = diagnostic
    )
}


## Generalized BCa (gbca) diagnostic.
## Tests adequacy of the BCa approximation by comparing with a warped
## normal model. If BCa is exact, the diagnostic D(z) is constant = 1.
## Efron (1987) Sec 7.
##
## Returns list with $equiv (2-row matrix: alpha, alphatil) and
## $Dmatrix (2-col matrix: zz, diag_smooth) and $sese (named vector).
gbca_diagnostic <- function(Y, tt, t0, a0, sdjack, local_dir, alpha) {
    B <- length(tt)
    m <- ncol(Y)
    sdboot <- stats::sd(tt)
    z0 <- stats::qnorm(sum(tt < t0) / B)

    ## Big-A raw acceleration via delta method. Efron (1987) Sec 7.
    std_tt <- (tt - mean(tt)) / sdboot
    A <- 0.5 * (sdboot / sdjack) * sum(std_tt^2 * local_dir) / B

    ## RMS covariance stability measure
    YY <- scale(Y, center = TRUE, scale = FALSE)
    Del <- stats::cov(YY, std_tt^2) / 2
    rms <- sqrt((m / (m - 1)) * sum(Del^2))
    sese <- c(A = A, rms = rms)

    ## Grid of standard normal quantiles for diagnostic evaluation
    zz <- seq(-2.5, 2.5, 0.05)
    cdf_grid <- stats::pnorm(zz)
    n_grid <- length(cdf_grid)
    tta <- as.vector(stats::quantile(tt, cdf_grid))

    ## ecdf_vals[i] = empirical CDF of local_dir up to bootstrap quantile
    ecdf_vals <- numeric(n_grid)
    for (i in seq_len(n_grid)) ecdf_vals[i] <- sum(local_dir[tt <= tta[i]]) / B

    ## H = ratio of empirical CDF to normal density; normalize at center
    H <- ecdf_vals / stats::dnorm(zz)
    i0 <- 51  # center index (z=0)
    diag_raw <- H / H[i0]
    diag_smooth <- stats::smooth.spline(zz, diag_raw, df = 5)$y

    Dmatrix <- cbind(zz, diag_smooth)

    ## Warping transformation
    scaled_a <- a0 / (1 - a0 * z0)
    I <- (cumsum(1 / diag_smooth) - 0.5 / diag_smooth) * diff(zz)[1]
    I <- I - I[i0]
    if (abs(scaled_a) < 1e-06)
        ww <- I else ww <- (exp(scaled_a * I) - 1) / scaled_a

    ## BCa transformation on original scale
    Z <- z0 + (z0 + zz) / (1 - a0 * (z0 + zz))
    bet <- stats::pnorm(Z)

    ## Warping functions
    wfun <- function(z) stats::approx(zz, ww, z, rule = 2, ties = mean)$y
    winv <- function(w) stats::approx(ww, zz, w, rule = 2, ties = mean)$y
    Phitil <- function(w) stats::pnorm(winv(w))

    ## Apply BCa formula on warped scale
    zt0 <- wfun(z0)
    zzt <- -wfun(rev(zz))
    Zt <- zt0 + (zt0 + zzt) / (1 - a0 * (zt0 + zzt))
    bett <- Phitil(Zt)

    ## equiv: alpha-tilde = gbca-adjusted coverage
    atil <- stats::approx(bett, cdf_grid, bet, rule = 2, ties = mean)$y
    altil <- stats::approx(cdf_grid, atil, alpha, rule = 2, ties = mean)$y
    equiv <- rbind(alpha, altil)
    dimnames(equiv)[[2]] <- rep(" ", length(alpha))
    dimnames(equiv)[[1]] <- c("alpha", "alphatil")

    list(equiv = equiv, Dmatrix = Dmatrix, sese = sese)
}
