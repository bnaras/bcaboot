#' Parametric BCa bootstrap confidence intervals
#'
#' Computes parametric bootstrap confidence intervals for a real-valued
#' parameter theta in a p-parameter exponential family.
#'
#' @param t0 Observed estimate of theta, usually by maximum likelihood.
#' @param tt A vector of B parametric bootstrap replications of theta.
#' @param bb A B by p matrix of natural sufficient vectors, where p is
#'   the dimension of the exponential family.
#' @param conf.level Confidence levels for the intervals.
#' @param n_jack Number of jackknife repetitions for internal SE.
#'   Set to 0 to skip.
#' @param jack_groups Number of groups per jackknife fold.
#' @param truncation Truncation parameter for acceleration calculation.
#'   Can be a vector to compute acceleration at multiple levels.
#' @param kl_fraction Proportion of "nearby" b vectors used in
#'   gradient estimation.
#' @param conf_density Logical; if TRUE, return BCa confidence density
#'   weights.
#' @param func Optional function for ABC (analytical bootstrap) limits.
#'
#' @return An object of class `"bcaboot"` with components:
#'   \describe{
#'     \item{limits}{9-row matrix of confidence limits}
#'     \item{stats}{2-row matrix of estimates and jackknife SEs}
#'     \item{B_mean}{Bootstrap sample size and mean of replicates}
#'     \item{ustats}{Bias-corrected estimator and its SE}
#'     \item{abc}{ABC limits and stats (when `func` provided)}
#'     \item{conf_density}{Confidence density weights (when requested)}
#'     \item{accel_matrix}{Acceleration at multiple truncation levels}
#'   }
#'
#' @references Efron B (1987). Better bootstrap confidence intervals.
#'   JASA 82, 171-200.
#' @references DiCiccio T and Efron B (1992). More accurate confidence
#'   intervals in exponential families. Biometrika, 231-245.
#' @references Efron B and Narasimhan B (2020). The Automatic
#'   Construction of Bootstrap Confidence Intervals. Journal of
#'   Computational and Graphical Statistics, 29(3), 608-619.
#'   \doi{10.1080/10618600.2020.1714633}
#'
#' @examples
#' data(diabetes, package = "bcaboot")
#' X <- diabetes$x
#' y <- scale(diabetes$y, center = TRUE, scale = FALSE)
#' lm.model <- lm(y ~ X - 1)
#' mu.hat <- lm.model$fitted.values
#' sigma.hat <- stats::sd(lm.model$residuals)
#' t0 <- summary(lm.model)$adj.r.squared
#' y.star <- sapply(mu.hat, rnorm, n = 500, sd = sigma.hat)
#' tt <- apply(y.star, 1, function(y) summary(lm(y ~ X - 1))$adj.r.squared)
#' b.star <- y.star %*% X
#' set.seed(1234)
#' bca_par(t0 = t0, tt = tt, bb = b.star)
#'
#' @export
bca_par <- function(t0, tt, bb,
                    conf.level = c(0.95, 0.90, 0.80, 0.68),
                    n_jack = 6,
                    jack_groups = 10,
                    truncation = 0.001,
                    kl_fraction = 0.333,
                    conf_density = FALSE,
                    func) {

    call <- match.call()

    ## Save rng state
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        stats::runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    alpha <- expand_alpha((1 - conf.level) / 2)
    B <- length(tt)
    J <- jack_groups
    K <- n_jack

    ## Core BCa computation via internal engine
    vl0 <- bca(t0, tt, bb = as.matrix(bb), alpha = alpha,
               trun = truncation, pct = kl_fraction)

    ## --- Remap bca() return to canonical structure ---
    ## bca() returns: $lims (rows: bcalims, %iles, Stand; cols: alpha),
    ##   $theta_a_z0_A_az, $sd_sdd_mean_B, $ustats
    theta_stats <- vl0$theta_a_z0_A_az  # c(theta, a, z0, A, az)
    sd_stats <- vl0$sd_sdd_mean_B       # c(sd, sdd, mean, B)

    ## Build canonical stats (subset matching the standard 5 columns)
    stats0 <- c(theta_stats["theta"], sd_stats["sd"],
                theta_stats["z0"], theta_stats["a"], NA_real_)
    names(stats0) <- c("theta", "sdboot", "z0", "a", "sdjack")

    ## Build canonical limits
    bca_lims <- vl0$lims[, "bcalims"]
    std_lims <- vl0$lims[, "Stand"]
    lims0 <- cbind(bca_lims, NA_real_, std_lims, NA_real_)
    dimnames(lims0) <- list(alpha, c("bca", "jacksd", "std", "pct"))

    ustats0 <- c(ustat = vl0$ustats[1], sdu = vl0$ustats[2])
    B_mean <- c(B, mean(tt))

    if (K == 0) {
        stats_mat <- rbind(stats0, jsd = rep(NA_real_, 5))
        dimnames(stats_mat) <- list(c("est", "jsd"),
                                    c("theta", "sdboot", "z0", "a", "sdjack"))
        result <- new_bcaboot(
            limits = lims0, stats = stats_mat, B_mean = B_mean,
            ustats = ustats0, call = call, seed = seed,
            method = "par", accel = NA_character_,
            accel_matrix = if (length(truncation) > 1) vl0$amat else NULL
        )
        if (!missing(func)) {
            vla <- abcpar(func, bb, alpha = alpha[alpha < 0.5])
            result$abc <- list(limits = vla$lims[, 1],
                               stats = vla$a.z0.cq[1:2])
        }
        return(result)
    }

    ## --- Internal SE via delete-d jackknife ---
    Limsd <- matrix(0, length(alpha), K)
    Thsd <- matrix(0, 5, K)
    Sdmsd <- matrix(0, 4, K)
    Ustm <- matrix(0, 2, K)

    for (k in seq_len(K)) {
        II <- matrix(sample(x = B, size = B), ncol = J)
        fold_lims <- matrix(0, length(alpha), J)
        fold_th <- matrix(0, 5, J)
        fold_sdm <- matrix(0, 4, J)
        fold_ust <- matrix(0, 2, J)
        for (j in seq_len(J)) {
            iij <- c(II[, -j])
            bbj <- as.matrix(bb[iij, ])
            ttj <- tt[iij]
            vlj <- bca(t0, ttj, bb = bbj, alpha = alpha,
                       trun = truncation, pct = kl_fraction)
            fold_lims[, j] <- vlj$lims[, 1]
            fold_th[, j] <- vlj$theta_a_z0_A_az
            fold_sdm[, j] <- vlj$sd_sdd_mean_B
            fold_ust[, j] <- vlj$ustats
        }
        Limsd[, k] <- apply(fold_lims, 1, stats::sd) * (J - 1) / sqrt(J)
        Thsd[, k] <- apply(fold_th, 1, stats::sd) * (J - 1) / sqrt(J)
        Sdmsd[, k] <- apply(fold_sdm, 1, stats::sd) * (J - 1) / sqrt(J)
        Ustm[, k] <- apply(fold_ust, 1, stats::sd) * (J - 1) / sqrt(J)
    }

    limsd <- rowMeans(Limsd)
    thsd <- rowMeans(Thsd)
    sdmsd <- rowMeans(Sdmsd)
    ustsd <- rowMeans(Ustm)

    ## Percentiles of bootstrap replicates at BCa limits
    pct <- numeric(length(alpha))
    for (i in seq_along(alpha)) pct[i] <- sum(tt <= lims0[i, "bca"]) / B

    ## Assemble final limits
    limits <- cbind(lims0[, "bca"], limsd, lims0[, "std"], pct)
    dimnames(limits) <- list(alpha, c("bca", "jacksd", "std", "pct"))

    ## Assemble final stats (canonical 5-column subset + jsd row)
    statsd <- c(sdmsd["sd"], sdmsd["sd"], thsd["z0"], thsd["a"], NA_real_)
    ## Actually, map correctly: theta_sd -> thsd[1], sdboot_sd -> sdmsd[1],
    ## z0_sd -> thsd[3], a_sd -> thsd[2], sdjack_sd -> NA
    statsd <- c(thsd[1], sdmsd[1], thsd[3], thsd[2], NA_real_)
    stats_mat <- rbind(stats0, statsd)
    dimnames(stats_mat) <- list(c("est", "jsd"),
                                c("theta", "sdboot", "z0", "a", "sdjack"))

    ## Ustats
    ustats_final <- c(ustat = ustats0["ustat"], sdu = ustats0["sdu"])

    ## Optional: accel_matrix for multi-truncation
    accel_matrix <- if (length(truncation) > 1) vl0$amat else NULL

    ## Optional: confidence density
    w <- NULL
    if (conf_density) {
        a_cd <- theta_stats["a"]
        z0_cd <- theta_stats["z0"]
        G <- (rank(tt) - 0.5) / B
        zth <- stats::qnorm(G) - z0_cd
        az <- 1 + a_cd * zth
        num <- stats::dnorm(zth / az - z0_cd)
        den <- az^2 * stats::dnorm(zth + z0_cd)
        w <- num / den
    }

    result <- new_bcaboot(
        limits = limits, stats = stats_mat, B_mean = B_mean,
        ustats = ustats_final, call = call, seed = seed,
        method = "par", accel = NA_character_,
        conf_density = w, accel_matrix = accel_matrix
    )

    ## Optional: ABC limits
    if (!missing(func)) {
        vla <- abcpar(func, bb, alpha = alpha[alpha < 0.5])
        result$abc <- list(limits = vla$lims[, 1],
                           stats = vla$a.z0.cq[1:2])
    }

    result
}
