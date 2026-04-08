#' Tidy a bcaboot result
#'
#' Returns a tibble with one row per (confidence level, method)
#' combination, following broom conventions.
#'
#' @param x A `bcaboot` object from [bca_nonpar()] or [bca_par()].
#' @param ... Ignored.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{conf.level}{Confidence level (e.g. 0.95)}
#'     \item{method}{`"bca"` or `"standard"` (or `"abc"` when available)}
#'     \item{estimate}{Point estimate (theta)}
#'     \item{conf.low}{Lower confidence limit}
#'     \item{conf.high}{Upper confidence limit}
#'     \item{jacksd.low}{Jackknife SE of lower limit (BCa only; `NA` otherwise)}
#'     \item{jacksd.high}{Jackknife SE of upper limit (BCa only; `NA` otherwise)}
#'   }
#'
#' @examples
#' data(diabetes, package = "bcaboot")
#' Xy <- cbind(diabetes$x, diabetes$y)
#' rfun <- function(Xy) {
#'   y <- Xy[, 11]; X <- Xy[, 1:10]
#'   summary(lm(y ~ X))$adj.r.squared
#' }
#' set.seed(1234)
#' result <- bca_nonpar(Xy, 1000, rfun, n_groups = 34, verbose = FALSE)
#' tidy(result)
#'
#' @export
tidy.bcaboot <- function(x, ...) {
    lims <- x$limits %||% x$lims
    stats <- x$stats

    ## Extract theta
    theta <- if (is.matrix(stats)) stats[1, 1] else stats[1]

    ## Read alpha from rownames, split into lower/upper
    alphas <- as.numeric(rownames(lims))
    lo_idx <- which(alphas < 0.5)
    hi_idx <- which(alphas > 0.5)
    conf_levels <- 1 - 2 * alphas[lo_idx]
    k <- length(conf_levels)

    ## Normalise column names for old-style objects
    cn <- colnames(lims)
    bca_col <- if ("bca" %in% cn) "bca" else cn[1]
    std_col <- if ("std" %in% cn) "std" else if ("Stand" %in% cn) "Stand" else NULL
    has_jacksd <- "jacksd" %in% cn && !all(is.na(lims[, "jacksd"]))

    ## BCa rows
    bca_df <- tibble::tibble(
        conf.level  = conf_levels,
        method      = "bca",
        estimate    = theta,
        conf.low    = lims[lo_idx, bca_col],
        conf.high   = lims[rev(hi_idx), bca_col],
        jacksd.low  = if (has_jacksd) lims[lo_idx, "jacksd"] else NA_real_,
        jacksd.high = if (has_jacksd) lims[rev(hi_idx), "jacksd"] else NA_real_
    )

    ## Standard rows
    if (!is.null(std_col)) {
        std_df <- tibble::tibble(
            conf.level  = conf_levels,
            method      = "standard",
            estimate    = theta,
            conf.low    = lims[lo_idx, std_col],
            conf.high   = lims[rev(hi_idx), std_col],
            jacksd.low  = NA_real_,
            jacksd.high = NA_real_
        )
        result <- rbind(bca_df, std_df)
    } else {
        result <- bca_df
    }

    ## ABC rows (from bca_par with func)
    abc <- x$abc
    if (!is.null(abc) && !is.null(abc$limits)) {
        abc_lims <- abc$limits
        abc_df <- tibble::tibble(
            conf.level  = conf_levels,
            method      = "abc",
            estimate    = theta,
            conf.low    = abc_lims[lo_idx],
            conf.high   = abc_lims[rev(hi_idx)],
            jacksd.low  = NA_real_,
            jacksd.high = NA_real_
        )
        result <- rbind(result, abc_df)
    }

    ## Sort by conf.level descending, then method
    result <- result[order(-result$conf.level, result$method), ]
    rownames(result) <- NULL
    result
}

#' Glance at a bcaboot result
#'
#' Returns a one-row tibble summarizing the bootstrap run.
#'
#' @param x A `bcaboot` object from [bca_nonpar()] or [bca_par()].
#' @param ... Ignored.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{method}{`"nonpar"` or `"par"`}
#'     \item{accel}{Acceleration method (`"regression"`, `"jackknife"`, or `NA`)}
#'     \item{theta}{Point estimate}
#'     \item{sdboot}{Bootstrap standard error}
#'     \item{z0}{Bias-correction constant}
#'     \item{a}{Acceleration constant}
#'     \item{sdjack}{Jackknife standard error (if available)}
#'     \item{B}{Number of bootstrap replicates}
#'     \item{boot_mean}{Mean of bootstrap replicates}
#'   }
#'
#' @examples
#' data(diabetes, package = "bcaboot")
#' Xy <- cbind(diabetes$x, diabetes$y)
#' rfun <- function(Xy) {
#'   y <- Xy[, 11]; X <- Xy[, 1:10]
#'   summary(lm(y ~ X))$adj.r.squared
#' }
#' set.seed(1234)
#' result <- bca_nonpar(Xy, 1000, rfun, n_groups = 34, verbose = FALSE)
#' glance(result)
#'
#' @export
glance.bcaboot <- function(x, ...) {
    stats <- x$stats
    B_mean <- x$B_mean %||% x$B.mean

    if (is.matrix(stats)) {
        row <- if ("est" %in% rownames(stats)) "est" else rownames(stats)[1]
        theta  <- stats[row, "theta"]
        sdboot <- if ("sdboot" %in% colnames(stats)) stats[row, "sdboot"] else
                  if ("sd" %in% colnames(stats)) stats[row, "sd"] else NA_real_
        z0     <- stats[row, "z0"]
        a      <- stats[row, "a"]
        sdjack <- if ("sdjack" %in% colnames(stats)) stats[row, "sdjack"] else NA_real_
    } else {
        theta  <- stats[1]
        sdboot <- stats[2]
        z0     <- stats["z0"]
        a      <- stats["a"]
        sdjack <- stats["sdjack"]
    }

    tibble::tibble(
        method    = x$method %||% NA_character_,
        accel     = x$accel %||% NA_character_,
        theta     = unname(theta),
        sdboot    = unname(sdboot),
        z0        = unname(z0),
        a         = unname(a),
        sdjack    = unname(sdjack),
        B         = B_mean[1],
        boot_mean = B_mean[2]
    )
}
