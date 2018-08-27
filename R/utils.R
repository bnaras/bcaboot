bcaboot.return <- function(x) {
    class(x) <- "bcaboot"
    x
}

#' @export
print.bcaboot <- function (x, digits = getOption("digits"), ...) {
    result <- x
    result$seed <- NULL
    print.default(result)
    invisible(x)
}
