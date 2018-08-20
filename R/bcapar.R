## Version of May 30, 2018
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
#' @param alpha Confidence values for which the confidence limits are
#'     calculated
#' @param J,K Parameters controlling the jackknife estimates of Monte
#'     Carlo error: `J` jackknife folds, with the jackknife standard
#'     errors averaged over `K` random divisions of `bb`
#' @param trun Truncation parameter used in the calculation of the
#'     acceleration `a`.
#' @param pct Proportion of "nearby" b vectors used in the calculation
#'     of `t.`, the gradient vector of theta.
#' @param cd If cd is 1 the bca confidence density is also returned
#'     (See Section 11.6 of "Computer Age Statistical Inference",
#'     B. Efron and T. Hastie.)
#' @param func Function `thetahat = t(b)`. If this is not missing then
#'     output includes _abc_ estimates (See
#'     "More accurate confidence intervals in exponential families",
#'     T. DiCiccio and B. Efron, Biometrika 1992 p231-245.)
#' @return a named list of several items:
#' \describe{
#' \item{`lims`}{bca confidence limits (first column) and the standard
#'     limits (fourth column). Also the abc limits (fifth column) if
#'     `func` is provided. The second column, `jacksd`, are the
#'     jackknife estimates of Monte Carlo error; `pctiles`, the third
#'     column are the proportion of the replicates `tt` less than each
#'     `bcalim` value}
#' \item{`stats`}{estimates and their jackknife Monte Carlo errors:
#'     `thet` = \eqn{\hat{\theta}}; `sd`, the bootstrap standard deviation
#'      for \eqn{\hat{\theta}}; `a` the acceleration estimate; `az` another
#'      acceleration estimate that depends less on extreme values of `tt`;
#'      `z0` the bias-correction estimate; `A` the big-A measure of raw
#'      acceleration; `sdd` delta method estimate for standard deviation of
#'      \eqn{\hat{\theta}}; `mean` the average of `tt`}
#' \item{abcstats}{The abc estimates of `a` and `z0`, returned if `func` was provided}
#' \item{ustats}{The bias-corrected estimator `2 * t0 - mean(tt)`. `ustats`
#'      gives `ustat`, an estimate `sdu` of its sampling error, and jackknife
#'      estimates of monte carlo error for both `ustat` and `sdu`. Also given
#'      is `B`, the number of bootstrap replications}
#' }
#' @references DiCiccio T and Efron B (1996). Bootstrap confidence
#'     intervals. Statistical Science 11, 189-228
#' @references T. DiCiccio and B. Efron. More accurate confidence intervals in exponential families.
#'     Biometrika (1992) p231-245.
#' @references Efron B (1987). Better bootstrap confidence intervals. JASA 82, 171-200
#' @references B. Efron and T. Hastie. Computer Age Statistical Inference. Cambridge University Press, 2016.
#'
#'
#' @import lars
#' @export
#' @examples
#' data(diabetes, package = "lars")
#' X <- diabetes$x
#' y <- scale(diabetes$y, center = TRUE, scale = FALSE)
#' lm.model <- lm(y ~ X - 1)
#' mu.hat <- lm.model$fitted.values
#' sigma.hat <- sd(lm.model$residuals)
#' t0 <- summary(lm.model)$adj.r.squared
#' y.star <- sapply(mu.hat, rnorm, n = 1000, sd = sigma.hat)
#' tt <- apply(y.star, 1, function(y) summary(lm(y ~ X - 1))$adj.r.squared)
#' b.star <- y.star %*% X
#' bcapar(t0 = t0, tt = tt, bb = b.star)
#' @export bcapar
bcapar <- function(t0, tt, bb, alpha = c(0.025, 0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95,
    0.975), J = 10, K = 6, trun = 0.001, pct = 0.333, cd = 0, func) {
    # if(!missing(func)) adds abc limits and stats K=0 skips jackknifing computes
    # confidence density weights 'w' if cd=1

    if (K == 0)
        return(bca(t0, tt, bb, alpha = alpha, trun = trun, pct = pct))
    call <- match.call()
    B <- length(tt)
    vl0 <- bca(t0, tt, bb = as.matrix(bb), alpha = alpha, trun = trun, pct = pct)
    Stand <- vl0$thet[1] + vl0$sd.m[1] * qnorm(alpha)
    Limsd <- matrix(0, length(alpha), K)
    Thsd <- matrix(0, 5, K)
    Sdmsd <- matrix(0, 4, K)
    Ustm <- matrix(0, 2, K)
    for (k in 1:K) {
        II <- sample(B, B)
        II <- matrix(II, ncol = J)
        lims <- matrix(0, length(alpha), J)
        th <- matrix(0, 5, J)
        sdm <- matrix(0, 4, J)
        usm <- matrix(0, 2, J)
        for (j in 1:J) {
            iij <- c(II[, -j])
            bbj <- as.matrix(bb[iij, ])
            ttj <- tt[iij]
            vlj <- bca(t0, ttj, bb = bbj, alpha = alpha, trun = trun, pct = pct)
            lims[, j] <- vlj$lims[, 1]
            th[, j] <- vlj$thet.a
            sdm[, j] <- vlj$sd.s
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
    dimnames(lim0) <- list(alpha, c("bcalims", "jacksd", "pctiles", "Stand"))
    th0 <- vl0$thet
    th0 <- rbind(th0, thsd)
    dimnames(th0) <- list(c("est", "jsd"), c("thet", "a", "z0", "A", "az"))
    sdm0 <- vl0$sd.s
    sdm0 <- rbind(sdm0, sdmsd)
    dimnames(sdm0) <- list(c("est", "jsd"), c("sd", "sdd", "mean", "B"))
    ust0 <- vl0$ustats
    ust0 <- rbind(ust0, ustsd)
    ust0 <- cbind(ust0, c(B, 0))
    dimnames(ust0) <- list(c("est", "jsd"), c("ustat", "sdu", "B"))
    ## lims <- round(lim0, roun)
    ## th0 <- round(th0, roun)
    ## sdm0 <- round(sdm0, roun)
    ## ust0 <- round(ust0, roun)
    lims <- lim0
    th0 <- th0
    stats <- cbind(th0, sdm0)
    stats <- stats[, c(1, 6, 2, 5, 3, 4, 7, 8)]

    vl <- list(call = call, lims = lims, stats = stats, ustats = ust0)
    if (length(trun) > 1)
        vl$amat <- vl0$amat
    if (cd == 1) {
        a <- stats[1, 4]
        z0 <- stats[1, 5]
        G <- (rank(tt) - 0.5)/B
        zth <- qnorm(G) - z0
        az <- 1 + a * zth
        num <- dnorm(zth/az - z0)
        den <- az^2 * dnorm(zth + z0)
        w <- num/den
        vl$w <- w
    }

    if (!missing(func)) {
        vla <- abcpar(func, bb, alpha = alpha[alpha < 0.5])
        abclims <- vla$lims[, 1]
        ## vl$lims <- round(cbind(vl$lims, abclims), roun)
        ## abcstats <- round(vla$a.z0.cq[1:2], roun)
        vl$lims <- cbind(vl$lims, abclims)
        abcstats <- vla$a.z0.cq[1:2]
        names(abcstats) <- c("a", "z0")
        vl$abcstats <- abcstats
    }
    vl
}
