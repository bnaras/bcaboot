## Shared internal helpers for BCa computation.

#' Compute BCa confidence limits from bootstrap replicates.
#'
#' BCa percentile formula: maps nominal alpha levels to bias-and-acceleration
#' adjusted percentiles of the bootstrap distribution. Efron (1987) Sec 2.
#'
#' @param z0 Bias-correction constant (scalar).
#' @param a Acceleration constant (scalar).
#' @param zalpha Normal quantiles at desired alpha levels (vector).
#' @param boot_reps Sorted vector of B bootstrap replicates.
#' @return Named list with `limits` (BCa quantiles), `bca_idx` (indices
#'   into sorted boot_reps), and `iles` (adjusted percentile levels).
#' @keywords internal
compute_bca_limits <- function(z0, a, zalpha, boot_reps) {
    B <- length(boot_reps)
    ## BCa percentile transformation. Efron (1987) Sec 2.
    iles <- stats::pnorm(z0 + (z0 + zalpha) / (1 - a * (z0 + zalpha)))
    bca_idx <- trunc(iles * B)
    bca_idx <- pmin(pmax(bca_idx, 1), B)
    sorted <- sort(boot_reps)
    list(limits = sorted[bca_idx], bca_idx = bca_idx, iles = iles)
}

#' Compute Poisson/multinomial deviance from uniform for bootstrap count matrix.
#'
#' Each row of Y is a bootstrap count vector. The deviance measures how far
#' each bootstrap sample departs from simple random sampling (uniform counts).
#' Limit when Yi=0: the term equals 2. EN20 appendix.
#'
#' @param Y B x m count matrix (rows = bootstrap samples, cols = observation counts).
#' @return Numeric vector of length B with deviance values.
#' @keywords internal
compute_kl_distance <- function(Y) {
    ## Vectorized: handle Yi=0 by replacing 0*log(0)=NaN with 0
    Y_safe <- Y
    Y_safe[Y == 0] <- 1  # log(1) = 0, so 2*1*log(1) = 0; then add back +2 below
    kl_term <- 2 * Y * log(Y_safe) - 2 * (Y - 1)
    ## Correct the Yi=0 entries: limit is 2
    kl_term[Y == 0] <- 2
    rowSums(kl_term)
}
