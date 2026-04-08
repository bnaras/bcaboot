# Plot BCa bootstrap confidence intervals

Creates a ggplot showing BCa and standard confidence intervals across
coverage levels.

## Usage

``` r
# S3 method for class 'bcaboot'
autoplot(object, ...)
```

## Arguments

- object:

  A `bcaboot` object from [`bca_nonpar()`](bca_nonpar.md) or
  [`bca_par()`](bca_par.md).

- ...:

  Ignored.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(diabetes, package = "bcaboot")
Xy <- cbind(diabetes$x, diabetes$y)
rfun <- function(Xy) {
  y <- Xy[, 11]; X <- Xy[, 1:10]
  summary(lm(y ~ X))$adj.r.squared
}
set.seed(1234)
result <- bca_nonpar(Xy, 1000, rfun, n_groups = 34, verbose = FALSE)
autoplot(result)
} # }
```
