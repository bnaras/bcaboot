## Version of May 30, 2018
##
## Parametric BCa confidence intervals for exponential families.
## Takes pre-computed bootstrap replications tt and sufficient statistic
## matrix bb (unlike bcajack/bcajack2 which do their own resampling).
## Delegates core BCa computation to internal bca() function.
## Optionally computes ABC (analytical) limits via abcpar().

#' Compute parametric bootstrap confidence intervals
#'
#' @description bcapar computes parametric bootstrap confidence
#'     intervals for a real-valued parameter theta in a p-parameter
#'     exponential family. It is described in Section 4 of the
#'     reference below.
#'
#' @inheritParams bcajack
#' @param t0 Observed estimate of theta, usually by maximum
#'     likelihood.
#' @param tt A vector of parametric bootstrap replications of theta of
#'     length `B`, usually large, say `B = 2000`
#' @param bb A `B` by `p` matrix of natural sufficient vectors, where
#'     `p` is the dimension of the exponential family.
#' @param J,K Parameters controlling the jackknife estimates of Monte
#'     Carlo error: `J` jackknife folds, with the jackknife standard
#'     errors averaged over `K` random divisions of `bb`
#' @param trun Truncation parameter used in the calculation of the
#'     acceleration `a`.
#' @param pct Proportion of "nearby" b vectors used in the calculation
#'     of `t.`, the gradient vector of theta.
#' @param cd If cd is 1 the bca confidence density is also returned;
#'     see Section 11.6 in reference Efron and Hastie (2016) below
#' @param func Function \eqn{\hat{\theta} = func(b)}. If this is not missing then
#'     output includes _abc_ estimates; see reference DiCiccio and Efron (1992) below
#' @return a named list of several items:
#'
#' * __lims__ : Bca confidence limits (first column) and the standard
#'     limits (fourth column). Also the abc limits (fifth column) if
#'     `func` is provided. The second column, `jacksd`, are the
#'     jackknife estimates of Monte Carlo error; `pct`, the third
#'     column are the proportion of the replicates `tt` less than each
#'     `bcalim` value
#'
#' * __stats__ : Estimates and their jackknife Monte Carlo errors:
#'     `theta` = \eqn{\hat{\theta}}; `sd`, the bootstrap standard deviation
#'      for \eqn{\hat{\theta}}; `a` the acceleration estimate; `az` another
#'      acceleration estimate that depends less on extreme values of `tt`;
#'      `z0` the bias-correction estimate; `A` the big-A measure of raw
#'      acceleration; `sdd` delta method estimate for standard deviation of
#'      \eqn{\hat{\theta}}; `mean` the average of `tt`
#'
#' * __abcstats__ : The abc estimates of `a` and `z0`, returned if `func` was provided
#'
#' * __ustats__ : The bias-corrected estimator `2 * t0 - mean(tt)`. `ustats`
#'      gives `ustat`, an estimate `sdu` of its sampling error, and jackknife
#'      estimates of monte carlo error for both `ustat` and `sdu`. Also given
#'      is `B`, the number of bootstrap replications
#'
#' * __seed__ : The random number state for reproducibility
#'
#' @references DiCiccio T and Efron B (1996). Bootstrap confidence
#'     intervals. Statistical Science 11, 189-228
#' @references T. DiCiccio and B. Efron. More accurate confidence intervals in exponential families.
#'     Biometrika (1992) p231-245.
#' @references Efron B (1987). Better bootstrap confidence intervals. JASA 82, 171-200
#' @references B. Efron and T. Hastie. Computer Age Statistical Inference. Cambridge University Press, 2016.
#' @references B. Efron and B. Narasimhan. Automatic Construction of Bootstrap Confidence Intervals, 2018.
#'
#' @export
#' @examples
#' data(diabetes, package = "bcaboot")
#' X <- diabetes$x
#' y <- scale(diabetes$y, center = TRUE, scale = FALSE)
#' lm.model <- lm(y ~ X - 1)
#' mu.hat <- lm.model$fitted.values
#' sigma.hat <- stats::sd(lm.model$residuals)
#' t0 <- summary(lm.model)$adj.r.squared
#' y.star <- sapply(mu.hat, rnorm, n = 1000, sd = sigma.hat)
#' tt <- apply(y.star, 1, function(y) summary(lm(y ~ X - 1))$adj.r.squared)
#' b.star <- y.star %*% X
#' set.seed(1234)
#' bcapar(t0 = t0, tt = tt, bb = b.star)
#' @export bcapar
bcapar <- function(t0, tt, bb,
                   alpha = c(0.025, 0.05, 0.1, 0.16),
                   J = 10, K = 6, trun = 0.001, pct = 0.333, cd = 0, func) {

    ## Save rng state
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        stats::runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    alpha <- expand_alpha(alpha)

    if (K == 0) {
        result <- bca(t0, tt, bb, alpha = alpha, trun = trun, pct = pct)
        result$seed <- seed
        bcaboot.return(result)
    }
    call <- match.call()
    B <- length(tt)
    vl0 <- bca(t0, tt, bb = as.matrix(bb), alpha = alpha, trun = trun, pct = pct)
    Stand <- vl0$theta_a_z0_A_az[1] + vl0$sd_sdd_mean_B[1] * stats::qnorm(alpha)
    Limsd <- matrix(0, length(alpha), K)
    Thsd <- matrix(0, 5, K)
    Sdmsd <- matrix(0, 4, K)
    Ustm <- matrix(0, 2, K)
    ## Internal SE via delete-d jackknife of the bootstrap.
    ## Split B replications into J groups, recompute bca() leaving one group
    ## out at a time. Repeat K times, average. Scale: (J-1)/sqrt(J).
    for (k in seq_len(K)) {
        II <- sample(x = B, size = B)
        II <- matrix(II, ncol = J)
        lims <- matrix(0, length(alpha), J)
        th <- matrix(0, 5, J)
        sdm <- matrix(0, 4, J)
        usm <- matrix(0, 2, J)
        for (j in seq_len(J)) {
            iij <- c(II[, -j])
            bbj <- as.matrix(bb[iij, ])
            ttj <- tt[iij]
            vlj <- bca(t0, ttj, bb = bbj, alpha = alpha, trun = trun, pct = pct)
            lims[, j] <- vlj$lims[, 1]
            th[, j] <- vlj$theta_a_z0_A_az
            sdm[, j] <- vlj$sd_sdd_mean_B
            usm[, j] <- vlj$ustats
        }

        Limsd[, k] <- apply(lims, 1, sd) * (J - 1)/sqrt(J)
        Thsd[, k] <- apply(th, 1, sd) * (J - 1)/sqrt(J)
        Sdmsd[, k] <- apply(sdm, 1, sd) * (J - 1)/sqrt(J)
        Ustm[, k] <- apply(usm, 1, sd) * (J - 1)/sqrt(J)
    }

    limsd <- rowMeans(Limsd)
    thsd <- rowMeans(Thsd)
    sdmsd <- rowMeans(Sdmsd)
    ustsd <- rowMeans(Ustm)
    lim0 <- vl0$lims
    lim0 <- cbind(lim0[, 1], limsd, lim0[, 2:3])
    ## Standardized column order: bca, jacksd, std, pct (matching bcajack/bcajack2)
    dimnames(lim0) <- list(alpha, c("bca", "jacksd", "pct", "std"))
    lim0 <- lim0[, c("bca", "jacksd", "std", "pct")]
    th0 <- vl0$theta_a_z0_A_az
    th0 <- rbind(th0, thsd)
    dimnames(th0) <- list(c("est", "jsd"), c("theta", "a", "z0", "A", "az"))
    sdm0 <- vl0$sd_sdd_mean_B
    sdm0 <- rbind(sdm0, sdmsd)
    dimnames(sdm0) <- list(c("est", "jsd"), c("sd", "sdd", "mean", "B"))
    ust0 <- vl0$ustats
    ust0 <- rbind(ust0, ustsd)
    ust0 <- cbind(ust0, c(B, 0))
    dimnames(ust0) <- list(c("est", "jsd"), c("ustat", "sdu", "B"))
    lims <- lim0
    stats <- cbind(th0, sdm0)
    stats <- stats[, c(1, 6, 2, 5, 3, 4, 7, 8)]

    vl <- list(call = call, lims = lims, stats = stats, ustats = ust0, seed = seed)
    if (length(trun) > 1)
        vl$amat <- vl0$amat
    ## Confidence density weights: ratio of BCa density to standard normal density.
    ## See Efron & Hastie (2016) Sec 11.6.
    if (cd == 1) {
        a <- stats[1, 4]
        z0 <- stats[1, 5]
        G <- (rank(tt) - 0.5)/B
        zth <- stats::qnorm(G) - z0
        az <- 1 + a * zth
        num <- stats::dnorm(zth/az - z0)
        den <- az^2 * stats::dnorm(zth + z0)
        w <- num/den
        vl$w <- w
    }

    if (!missing(func)) {
        vla <- abcpar(func, bb, alpha = alpha[alpha < 0.5])
        abc <- vla$lims[, 1]
        vl$lims <- cbind(vl$lims, abc)
        abcstats <- vla$a.z0.cq[1:2]
        names(abcstats) <- c("a", "z0")
        vl$abcstats <- abcstats
    }
    bcaboot.return(vl)
}
