# Deprecated functions --------------------------------------------------------
#
# These wrappers provide backward compatibility with the pre-0.4 API.
# They translate old parameter names to the new API and emit a
# once-per-session deprecation warning via lifecycle.

#' @title Nonparametric BCa (deprecated)
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `bcajack()` was deprecated in bcaboot 1.0 in favour of [bca_nonpar()].
#' @inheritParams bca_nonpar
#' @param m Number of groups (now `n_groups`).
#' @param mr Group repetitions (now `group_reps`).
#' @param K Jackknife repetitions (now `n_jack`).
#' @param J Jackknife groups (now `jack_groups`).
#' @param alpha Alpha levels (now use `conf.level`).
#' @keywords internal
#' @export
bcajack <- function(x, B, func, ..., m = nrow(x), mr = 5, K = 2, J = 10,
                    alpha = c(0.025, 0.05, 0.1, 0.16), verbose = TRUE) {
    lifecycle::deprecate_warn("1.0", "bcajack()", "bca_nonpar()")
    ## Handle B as pre-computed vector (old bcajack overloading)
    boot_data <- NULL
    if (length(B) > 1) {
        boot_data <- list(tt = B, t0 = func(x, ...))
        B <- length(boot_data$tt)
    }
    bca_nonpar(x, B, func, ...,
               accel = "jackknife",
               conf.level = 1 - 2 * alpha[alpha < 0.5],
               n_groups = m, group_reps = mr,
               n_jack = K, jack_groups = J,
               boot_data = boot_data,
               verbose = verbose)
}

#' @title Nonparametric BCa via regression (deprecated)
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `bcajack2()` was deprecated in bcaboot 1.0 in favour of [bca_nonpar()].
#' @inheritParams bca_nonpar
#' @param m Number of groups (now `n_groups`).
#' @param mr Unused (included for compatibility).
#' @param pct KL fraction (now `kl_fraction`).
#' @param K Jackknife repetitions (now `n_jack`).
#' @param J Jackknife groups (now `jack_groups`).
#' @param alpha Alpha levels (now use `conf.level`).
#' @keywords internal
#' @export
bcajack2 <- function(x, B, func, ..., m = NULL, mr, pct = 0.333, K = 2, J = 12,
                     alpha = c(0.025, 0.05, 0.1, 0.16), verbose = TRUE) {
    lifecycle::deprecate_warn("1.0", "bcajack2()", "bca_nonpar()")
    if (!missing(mr))
        warning("'mr' argument is unused in bcajack2/bca_nonpar with regression acceleration")
    ## Handle B as Blist
    boot_data <- if (is.list(B)) B else NULL
    if (is.null(m) && !is.list(B)) m <- nrow(x)
    bca_nonpar(x, if (is.list(B)) length(B$tt) else B, func, ...,
               accel = "regression",
               conf.level = 1 - 2 * alpha[alpha < 0.5],
               n_groups = m, kl_fraction = pct,
               n_jack = K, jack_groups = J,
               boot_data = boot_data,
               verbose = verbose)
}

#' @title Parametric BCa (deprecated)
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `bcapar()` was deprecated in bcaboot 1.0 in favour of [bca_par()].
#' @inheritParams bca_par
#' @param alpha Alpha levels (now use `conf.level`).
#' @param K Jackknife repetitions (now `n_jack`).
#' @param J Jackknife groups (now `jack_groups`).
#' @param trun Truncation (now `truncation`).
#' @param pct KL fraction (now `kl_fraction`).
#' @param cd Confidence density flag (now `conf_density`).
#' @keywords internal
#' @export
bcapar <- function(t0, tt, bb,
                   alpha = c(0.025, 0.05, 0.1, 0.16),
                   J = 10, K = 6, trun = 0.001, pct = 0.333, cd = 0, func) {
    lifecycle::deprecate_warn("1.0", "bcapar()", "bca_par()")
    if (missing(func)) {
        bca_par(t0, tt, bb,
                conf.level = 1 - 2 * alpha[alpha < 0.5],
                n_jack = K, jack_groups = J,
                truncation = trun, kl_fraction = pct,
                conf_density = as.logical(cd))
    } else {
        bca_par(t0, tt, bb,
                conf.level = 1 - 2 * alpha[alpha < 0.5],
                n_jack = K, jack_groups = J,
                truncation = trun, kl_fraction = pct,
                conf_density = as.logical(cd),
                func = func)
    }
}

#' @title BCa plot (deprecated)
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `bcaplot()` was deprecated in bcaboot 1.0 in favour of
#' `autoplot()`.
#' @param vl Output of a bcaboot computation function.
#' @param ... Additional arguments (ignored).
#' @keywords internal
#' @export
bcaplot <- function(vl, ...) {
    lifecycle::deprecate_warn("1.0", "bcaplot()", "autoplot()")
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        print(ggplot2::autoplot(vl, ...))
    } else {
        warning("ggplot2 is required for plotting. ",
                "Install it with install.packages('ggplot2').",
                call. = FALSE)
    }
    invisible(vl)
}
