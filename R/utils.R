## Canonical constructor for bcaboot objects.
##
## Guarantees a consistent structure regardless of which function produced it.
## All computation functions should call this to build the return object.
##
## @param limits 9-row x 4-column matrix (alpha x c("bca","jacksd","std","pct")).
##   NA_real_ for columns not computed.
## @param stats 2-row x 5-column matrix (c("est","jsd") x c("theta","sdboot","z0","a","sdjack")).
##   "jsd" row is NA_real_ when n_jack = 0.
## @param B_mean length-2 numeric c(B, mean(tt)).
## @param ustats named vector c(ustat, sdu), or NA when not computable.
## @param call the matched call.
## @param seed .Random.seed at entry.
## @param method "nonpar" or "par".
## @param accel "regression", "jackknife", or NA for parametric.
## @param diagnostic list with $equiv, $Dmatrix (gbca); NULL when absent.
## @param abc list with $limits, $stats (abc method); NULL when absent.
## @param conf_density numeric vector of confidence density weights; NULL when absent.
## @param accel_matrix matrix of acceleration at multiple truncation levels; NULL when absent.
## @param sese named vector c(A, rms); NULL when absent.
new_bcaboot <- function(limits, stats, B_mean, ustats, call, seed,
                        method = "nonpar", accel = "regression",
                        diagnostic = NULL, abc = NULL,
                        conf_density = NULL, accel_matrix = NULL,
                        sese = NULL) {
    structure(
        list(
            limits        = limits,
            stats         = stats,
            B_mean        = B_mean,
            ustats        = ustats,
            method        = method,
            accel         = accel,
            call          = call,
            seed          = seed,
            diagnostic    = diagnostic,
            abc           = abc,
            conf_density  = conf_density,
            accel_matrix  = accel_matrix,
            sese          = sese
        ),
        class = "bcaboot"
    )
}

## Legacy constructor — used by existing functions until they are rewritten in Phase 3.
bcaboot.return <- function(x) {
    class(x) <- "bcaboot"
    x
}

#' @export
`$.bcaboot` <- function(x, name) {
    aliases <- c(lims = "limits", B.mean = "B_mean")
    if (name %in% names(aliases)) {
        .subset2(x, aliases[[name]]) %||% .subset2(x, name)
    } else {
        .subset2(x, name)
    }
}

#' @export
print.bcaboot <- function(x, digits = getOption("digits"), ...) {
    ## Detect whether this is a new-style (has $limits) or old-style (has $lims) object
    lims <- x$limits %||% x$lims
    stats <- x$stats

    ## --- Header ---
    method_str <- x$method
    accel_str  <- x$accel
    cli::cli_h3("BCa Bootstrap Confidence Intervals")
    if (!is.null(method_str)) {
        if (!is.null(accel_str) && !is.na(accel_str))
            cli::cli_text("Method: {method_str} ({accel_str} acceleration)")
        else
            cli::cli_text("Method: {method_str}")
    }

    ## Extract theta and sdboot from stats (vector or matrix)
    if (is.matrix(stats)) {
        theta  <- stats[1, "theta"]
        sdboot <- stats[1, grep("^sd", colnames(stats))[1]]
    } else {
        theta  <- stats[1]
        sdboot <- stats[2]
    }
    B_mean <- x$B_mean %||% x$B.mean
    cli::cli_text("B = {B_mean[1]}, theta = {format(theta, digits = digits)}, sdboot = {format(sdboot, digits = digits)}")
    cli::cli_text("")

    ## --- Confidence limits table ---
    ## Read alpha from rownames
    alphas <- as.numeric(rownames(lims))
    lo_idx <- which(alphas < 0.5)
    hi_idx <- which(alphas > 0.5)

    if (length(lo_idx) > 0 && length(lo_idx) == length(hi_idx)) {
        conf_levels <- 1 - 2 * alphas[lo_idx]
        bca_col <- if ("bca" %in% colnames(lims)) "bca" else colnames(lims)[1]
        std_col <- if ("std" %in% colnames(lims)) "std" else if ("Stand" %in% colnames(lims)) "Stand" else NULL

        tbl <- data.frame(
            conf.level = conf_levels,
            bca.lo     = lims[lo_idx, bca_col],
            bca.hi     = lims[rev(hi_idx), bca_col]
        )
        if (!is.null(std_col)) {
            tbl$std.lo <- lims[lo_idx, std_col]
            tbl$std.hi <- lims[rev(hi_idx), std_col]
        }

        cli::cli_text("{.strong Confidence limits:}")
        print(tbl, digits = digits, row.names = FALSE)
    } else {
        cli::cli_text("{.strong Limits:}")
        print(lims, digits = digits)
    }

    ## --- Diagnostics line ---
    cli::cli_text("")
    cli::cli_text("{.strong Diagnostics:}")
    if (is.matrix(stats)) {
        row <- if ("est" %in% rownames(stats)) "est" else rownames(stats)[1]
        z0 <- stats[row, "z0"]
        a  <- stats[row, "a"]
        sdjack <- if ("sdjack" %in% colnames(stats)) stats[row, "sdjack"] else NA
    } else {
        z0     <- stats["z0"]
        a      <- stats["a"]
        sdjack <- stats["sdjack"]
    }
    diag_parts <- c(
        paste0("z0 = ", format(z0, digits = digits)),
        paste0("a = ", format(a, digits = digits))
    )
    if (!is.na(sdjack))
        diag_parts <- c(diag_parts, paste0("sdjack = ", format(sdjack, digits = digits)))
    cli::cli_text("  {paste(diag_parts, collapse = ', ')}")

    invisible(x)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
