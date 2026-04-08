#' Plot BCa bootstrap confidence intervals
#'
#' Creates a ggplot showing BCa and standard confidence intervals
#' across coverage levels.
#'
#' @param object A `bcaboot` object from [bca_nonpar()] or [bca_par()].
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' data(diabetes, package = "bcaboot")
#' Xy <- cbind(diabetes$x, diabetes$y)
#' rfun <- function(Xy) {
#'   y <- Xy[, 11]; X <- Xy[, 1:10]
#'   summary(lm(y ~ X))$adj.r.squared
#' }
#' set.seed(1234)
#' result <- bca_nonpar(Xy, 1000, rfun, n_groups = 34, verbose = FALSE)
#' autoplot(result)
#' }
#'
#' @export
autoplot.bcaboot <- function(object, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        cli::cli_abort(c(
            "Package {.pkg ggplot2} is required for {.fun autoplot}.",
            "i" = "Install it with {.code install.packages(\"ggplot2\")}."
        ))
    }

    df <- tidy.bcaboot(object)
    theta <- df$estimate[1]

    ## Split into BCa and other methods for different visual treatment
    bca_df <- df[df$method == "bca", ]
    other_df <- df[df$method != "bca", ]

    ## Offset for non-BCa methods so segments don't overlap
    dx <- if (nrow(bca_df) > 1) diff(range(bca_df$conf.level)) * 0.025 else 0.02

    ## Coverage labels for legend
    cov_labels <- paste0(round(100 * bca_df$conf.level), "%")
    bca_df$coverage <- factor(cov_labels, levels = cov_labels)

    p <- ggplot2::ggplot() +

        ## BCa vertical segments colored by coverage
        ggplot2::geom_segment(
            data = bca_df,
            ggplot2::aes(x = .data$conf.level, xend = .data$conf.level,
                         y = .data$conf.low, yend = .data$conf.high,
                         colour = .data$coverage),
            linewidth = 1.5, lineend = "round"
        ) +

        ## BCa envelope lines
        ggplot2::geom_line(
            data = bca_df,
            ggplot2::aes(x = .data$conf.level, y = .data$conf.low),
            colour = "grey30", linewidth = 0.4
        ) +
        ggplot2::geom_line(
            data = bca_df,
            ggplot2::aes(x = .data$conf.level, y = .data$conf.high),
            colour = "grey30", linewidth = 0.4
        ) +

        ## BCa endpoint points
        ggplot2::geom_point(
            data = bca_df,
            ggplot2::aes(x = .data$conf.level, y = .data$conf.low,
                         colour = .data$coverage),
            size = 3
        ) +
        ggplot2::geom_point(
            data = bca_df,
            ggplot2::aes(x = .data$conf.level, y = .data$conf.high,
                         colour = .data$coverage),
            size = 3
        ) +

        ## Theta reference line
        ggplot2::geom_hline(yintercept = theta, linewidth = 0.5,
                            linetype = "dotted", colour = "#8E44AD") +
        ggplot2::annotate(
            "text", x = max(bca_df$conf.level), y = theta,
            label = paste0("hat(theta) == ", signif(theta, 4)),
            parse = TRUE, hjust = 1, vjust = -0.5,
            size = 3.3, colour = "#8E44AD"
        ) +

        ## Scales
        ggplot2::scale_x_continuous(
            breaks = bca_df$conf.level,
            labels = cov_labels,
            expand = ggplot2::expansion(mult = c(0.05, 0.12))
        ) +
        ggplot2::labs(
            x = "Coverage", y = "Confidence limit", colour = "Coverage",
            caption = "Solid colored = BCa \u2022 Grey dashed = Standard"
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(
            legend.position    = "right",
            panel.grid.minor   = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            plot.caption       = ggplot2::element_text(colour = "grey50",
                                                       size = 9, hjust = 0),
            legend.title       = ggplot2::element_text(face = "bold", size = 10)
        )

    ## Standard / ABC method lines (offset, grey dashed)
    if (nrow(other_df) > 0) {
        for (meth in unique(other_df$method)) {
            mdf <- other_df[other_df$method == meth, ]
            lt <- if (meth == "standard") "dashed" else "dotdash"

            p <- p +
                ggplot2::geom_segment(
                    data = mdf,
                    ggplot2::aes(x = .data$conf.level + dx,
                                 xend = .data$conf.level + dx,
                                 y = .data$conf.low, yend = .data$conf.high),
                    linewidth = 0.7, linetype = lt, colour = "grey50"
                ) +
                ggplot2::geom_line(
                    data = mdf,
                    ggplot2::aes(x = .data$conf.level + dx, y = .data$conf.low),
                    colour = "grey50", linewidth = 0.4, linetype = lt
                ) +
                ggplot2::geom_line(
                    data = mdf,
                    ggplot2::aes(x = .data$conf.level + dx, y = .data$conf.high),
                    colour = "grey50", linewidth = 0.4, linetype = lt
                ) +
                ggplot2::geom_point(
                    data = mdf,
                    ggplot2::aes(x = .data$conf.level + dx, y = .data$conf.low),
                    colour = "grey50", size = 1.8, shape = 17
                ) +
                ggplot2::geom_point(
                    data = mdf,
                    ggplot2::aes(x = .data$conf.level + dx, y = .data$conf.high),
                    colour = "grey50", size = 1.8, shape = 17
                )
        }
    }

    ## Jackknife SE error bars (BCa only, when available)
    has_jacksd <- any(!is.na(bca_df$jacksd.low))
    if (has_jacksd) {
        jdf_lo <- bca_df[!is.na(bca_df$jacksd.low), ]
        jdf_hi <- bca_df[!is.na(bca_df$jacksd.high), ]

        p <- p +
            ggplot2::geom_errorbar(
                data = jdf_lo,
                ggplot2::aes(x = .data$conf.level,
                             ymin = .data$conf.low - .data$jacksd.low,
                             ymax = .data$conf.low + .data$jacksd.low),
                width = dx * 0.8, colour = "#C0392B", linewidth = 0.6
            ) +
            ggplot2::geom_errorbar(
                data = jdf_hi,
                ggplot2::aes(x = .data$conf.level,
                             ymin = .data$conf.high - .data$jacksd.high,
                             ymax = .data$conf.high + .data$jacksd.high),
                width = dx * 0.8, colour = "#C0392B", linewidth = 0.6
            )
    }

    p
}
