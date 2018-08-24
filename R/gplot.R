#' @title Plots of bca confidence limits
#'
#' @description \code{ggbcaplot} uses the output of \code{bcajack},
#'     \code{bcajack2}, or \code{bcapar} to plot bca and standard
#'     confidence limits for the parameter of interest.
#'
#' @param vl output of \code{bcajack}, \code{bcajack2}, or
#'     \code{bcapar}
#' @param main the main caption (can be empty)
#' @param xlab x axis label (supplied if not specified)
#' @param ylab y axis labels (supplied if not specified)
#' @param two-sided coverages are \eqn{1-2\alpha},
#'     e.g. \code{alpha=c(.025,.05)} plots intervals
#'     \code{[.025,.975]} and \code{[.05,.95]}. Default is
#'     \code{alpha=c(.025,.05,.1,.16)} giving coverages
#'     .95,.90,.80,.68
#' @param ... further args for ggplot
#'
#' @details confidence interval endpoints are plotted vertically versus
#'     two-sided coverages \eqn{1-2\alpha}. Bca limits in black,
#'     Standard limits in green (dashed.). If \code{vl$lims} includes
#'     the column "jacksd" of jackknife internal standard deviations
#'     then these are indicated by vertical red bars centered at the
#'     bca limit points.
#'
#' @importFrom ggplot2 geom_line geom_segment geom_point aes labs
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter
#' @importFrom tidyr gather spread
#' @export
ggbcaplot <- function(vl, main = "", xlab = "coverage", ylab = "limits",
                      alpha = c(0.025, 0.05, 0.1, 0.16), ...) {
    n_alpha <- length(alpha)
    alp <- rev(1 - 2 * alpha)
    vm <- vl$lims
    thet <- if (is.null(vl$thet)) vl$stat[1, 1] else vl$thet[1]

    vm <- data.frame(cbind(vm, alpha = as.numeric(rownames(vm))), stringsAsFactors = FALSE)
    vm <- vm[c(match(alpha, vm$alpha), match(1 - alpha, vm$alpha)), ]
    vm$type <- c(rep("lower", n_alpha), rep("upper", n_alpha))
    vm$alpha <- rep(vm$alpha[seq_len(n_alpha)], 2)
    vm$coverage <- 1 - 2 * vm$alpha
    vm %>%
        dplyr::select(-alpha) %>%
        tidyr::gather(key = key, value = value, -coverage, -type) %>%
        tidyr::spread(key = key, value = value) ->
        vmf
    lower <- vmf %>% dplyr::filter(type == "lower")
    upper <- vmf %>% dplyr::filter(type == "upper")

    g <- ggplot2::ggplot(...) +
        ggplot2::geom_line(mapping = ggplot2::aes(x = coverage, y = bcalims), data = lower) +
        ggplot2::geom_line(mapping = ggplot2::aes(x = coverage, y = bcalims), data = upper) +
        ggplot2::geom_point(mapping = ggplot2::aes(x = coverage, y = bcalims), data = lower) +
        ggplot2::geom_point(mapping = ggplot2::aes(x = coverage, y = bcalims), data = upper) +
        ggplot2::geom_line(mapping = ggplot2::aes(x = coverage, y = standard), data = lower, color = "blue", linetype = 3) +
        ggplot2::geom_line(mapping = ggplot2::aes(x = coverage, y = standard), data = upper, color = "blue", linetype = 3) +
        ggplot2::geom_point(mapping = ggplot2::aes(x = coverage, y = standard), data = lower, color = "blue") +
        ggplot2::geom_point(mapping = ggplot2::aes(x = coverage, y = standard), data = upper, color = "blue")

    g <- g +
        ggplot2::geom_segment(mapping = ggplot2::aes(x = coverage, xend = coverage, y = bcalims - jacksd, yend = bcalims + jacksd),
                     color = "red", data = vmf)
    g <- g +
        ggplot2::geom_hline(yintercept = thet, linetype = "dashed", color = "red")

    g <-  g +
        ggplot2::geom_segment(mapping = ggplot2::aes(x = coverage, xend = coverage, y = bcalims[2], yend = bcalims[1], group = coverage),
                              color = "green", data = vmf)

    g +
        labs(x = xlab, y = ylab, main = main)
}
