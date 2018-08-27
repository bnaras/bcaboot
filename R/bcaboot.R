#' Automatic Construction of Bootstrap Confidence Intervals
#'
#' Bootstrap confidence intervals depend on three elements: (a) the
#' cumulative distribution of the bootstrap replications, (b) the
#' bias-correction, and (c) the acceleration number that measures the
#' rate of change in the standard deviation of the estimate as the
#' data changes.  The first two of these depend only on the bootstrap
#' distribution, and not how it is generated: parametrically or
#' non-parametrically. Therefore, the only difference in a parametric
#' bca analysis would lie in the nonparametric estimation of the
#' acceleration, often a negligible error.
#'
#' @docType package
#' @name bcaboot
#'
NULL


#' Blood and other measurements in diabetics
#'
#' The \code{diabetes} data frame has 442 rows and 3 columns. These are the
#' data used in the Efron et al "Least Angle Regression" paper.
#'
#' The x matrix has been standardized to have unit L2 norm in each column and
#' zero mean. The matrix x2 consists of x plus certain interactions.
#'
#' @format This data frame contains the following columns:
#' - __x__ a matrix with 10 columns
#' - __y__ a numeric vector
#' - __x2__ a matrix with 64 columns
#'
#' @references Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle
#' Regression" (with discussion) \emph{Annals of Statistics}
#' @source
#' \url{https://web.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf}
#' @keywords datasets
#' @docType data
#' @name diabetes
NULL
