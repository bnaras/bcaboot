## Shared internal helpers for BCa computation.

## Expand user-supplied alpha values (lower tail only) into the full
## symmetric set including the median. E.g. c(0.025, 0.05, 0.1, 0.16) ->
## c(0.025, 0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95, 0.975)
expand_alpha <- function(alpha) {
    alpha <- alpha[alpha < 0.5]
    c(alpha, 0.5, rev(1 - alpha))
}

## Bias-corrected point estimate and its SE.
## ustat = 2*t0 - mean(tt) is the bootstrap bias-corrected estimator.
## sdu = SE from variance-stabilized influence: 2*primary - covariance-based.
## adj_infl must be pre-computed as (2 * primary_infl - cov_infl).
compute_ustats <- function(t0, tt, adj_infl, m) {
    ustat <- 2 * t0 - mean(tt)
    sdu <- sqrt(sum(adj_infl^2)) / m
    c(ustat = ustat, sdu = sdu)
}

## Regression-based acceleration estimation from bootstrap count matrix.
## Used by bcajack2 and bcanon (inside their qbca2 inner functions).
## Estimates influence via local regression on count vectors nearest
## to uniform (1,...,1). EN20 appendix.
##
## Returns list(a, sdjack, reg_infl, local_dir).
## local_dir is the normalized projection direction (needed by bcanon's gbca diagnostic).
regression_accel <- function(Y, tt, t0, pct) {
    m <- ncol(Y)
    B <- nrow(Y)

    ## KL distance: select pct fraction of bootstrap samples closest to uniform
    kl_dist <- compute_kl_distance(Y)
    kl_cutoff <- stats::quantile(kl_dist, pct)
    nearby_idx <- seq_len(B)[kl_dist <= kl_cutoff]

    if (length(nearby_idx) < m) {
        stop(sprintf(
            paste0("Regression for acceleration estimation is underdetermined: ",
                   "%d nearby bootstrap samples but %d columns in the count matrix. ",
                   "Possible fixes: (1) increase B so that B * kl_fraction > %d, ",
                   "(2) increase kl_fraction (currently %.3f), or ",
                   "(3) use accel = \"jackknife\" which avoids this regression entirely."),
            length(nearby_idx), m, m, pct
        ), call. = FALSE)
    }

    ## Influence via local regression on nearby count vectors
    reg_infl <- as.vector(m * stats::lm(tt[nearby_idx] ~ Y[nearby_idx, ] - 1)$coef)
    reg_infl <- reg_infl - mean(reg_infl)

    ## Acceleration: skewness of influence distribution. Efron (1987) Sec 6
    a <- (1/6) * sum(reg_infl^3) / sum(reg_infl^2)^1.5

    ## Jackknife SE from regression-based influences
    sdjack <- sqrt(sum(reg_infl^2)) / (m - 1)

    ## local_dir: normalized projection of centered counts onto influence direction
    ## Needed for bcanon's gbca diagnostic. Efron (1987) Sec 7.
    YY <- scale(Y, center = TRUE, scale = FALSE)
    local_dir <- as.vector(YY %*% reg_infl) / (m - 1)
    local_dir <- local_dir / stats::sd(local_dir)

    list(a = a, sdjack = sdjack, reg_infl = reg_infl, local_dir = local_dir)
}

## Classical delete-one (or delete-group) jackknife for acceleration estimation.
## Used by bcajack. Returns list(a, sdjack, jack_infl).
##
## When m == n: delete-one jackknife.
## When m < n: grouped jackknife, averaged over mr random partitions.
jackknife_accel <- function(x, func, ..., m, mr) {
    n <- nrow(x)
    u <- numeric(length = m)
    jk_norm <- sqrt(m * (m - 1))

    if (m == n) {
        for (i in seq_len(n)) {
            u[i] <- func(x[-i, , drop = FALSE], ...)
        }
        jack_infl <- (mean(u) - u) * (m - 1)
        a <- (1 / 6) * sum(jack_infl^3) / (sum(jack_infl^2))^1.5
        sdjack <- sqrt(sum(jack_infl^2)) / jk_norm
    } else if (m < n) {
        aa <- ssj <- numeric(mr)
        r <- n %% m
        seq_len_m <- seq_len(m)
        for (k in seq_len(mr)) {
            Imat <- matrix(sample.int(n, size = n - r), nrow = m)
            Iout <- setdiff(seq_len(n), Imat)
            for (j in seq_len_m) {
                Ij <- setdiff(seq_len_m, j)
                ij <- c(c(Imat[Ij, ], Iout))
                u[j] <- func(x[ij, , drop = FALSE])
            }
            jack_infl <- (mean(u) - u) * (m - 1)
            aa[k] <- (1/6) * sum(jack_infl^3)/(sum(jack_infl^2))^1.5
            ssj[k] <- sqrt(sum(jack_infl^2))/jk_norm
        }
        a <- mean(aa)
        sdjack <- mean(ssj)
        ## Use the last partition's influence values (consistent with original bcajack)
        jack_infl <- (mean(u) - u) * (m - 1)
    } else {
        stop("m must be <= n")
    }

    list(a = a, sdjack = sdjack, jack_infl = jack_infl)
}

## Bootstrap resampling with count matrix construction.
## Used by bcajack2 and bcanon. Returns list(tt, Y, t0).
##
## When boot_data is a list with components Y, tt, t0: uses those directly.
## When boot_data is NULL: generates B bootstrap resamples from x using func.
## Handles both m == n (standard) and m < n (grouped) resampling.
bootstrap_resample <- function(x, B, func, ..., m, verbose, boot_data = NULL) {
    if (!is.null(boot_data)) {
        return(list(Y = boot_data$Y, tt = boot_data$tt, t0 = boot_data$t0,
                    B = length(boot_data$tt)))
    }

    if (is.vector(x))
        x <- as.matrix(x)
    n <- nrow(x)
    if (is.null(m)) m <- n
    tt <- numeric(B)
    t0 <- func(x, ...)

    if (m == n) {
        ii <- sample(x = seq_len(n), size = n * B, replace = TRUE)
        ii <- matrix(ii, B)
        Y <- matrix(0, B, n)
        if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
        for (k in seq_len(B)) {
            ik <- ii[k, ]
            tt[k] <- func(x[ik, , drop = FALSE], ...)
            Y[k, ] <- table(c(ik, 1:n)) - 1
            if (verbose) utils::setTxtProgressBar(pb, k)
        }
        if (verbose) close(pb)
    } else if (m < n) {
        r <- n %% m
        Imat <- matrix(sample(x = seq_len(n), size = n - r), m)
        Iout <- setdiff(seq_len(n), Imat)
        ii <- sample(x = seq_len(m), size = m * B, replace = TRUE)
        ii <- matrix(ii, B)
        Y <- matrix(0, B, m)
        if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
        for (k in seq_len(B)) {
            ik <- ii[k, ]
            Ik <- c(t(Imat[ik, ]))
            Ik <- c(Ik, Iout)
            tt[k] <- func(x[Ik, , drop = FALSE], ...)
            Y[k, ] <- table(c(ik, 1:m)) - 1
            if (verbose) utils::setTxtProgressBar(pb, k)
        }
        if (verbose) close(pb)
    } else {
        stop("m must be <= n")
    }

    list(Y = Y, tt = tt, t0 = t0, B = B)
}

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
