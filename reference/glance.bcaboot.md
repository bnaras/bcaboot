# Glance at a bcaboot result

Returns a one-row tibble summarizing the bootstrap run.

## Usage

``` r
# S3 method for class 'bcaboot'
glance(x, ...)
```

## Arguments

- x:

  A `bcaboot` object from [`bca_nonpar()`](bca_nonpar.md) or
  [`bca_par()`](bca_par.md).

- ...:

  Ignored.

## Value

A tibble with columns:

- method:

  `"nonpar"` or `"par"`

- accel:

  Acceleration method (`"regression"`, `"jackknife"`, or `NA`)

- theta:

  Point estimate

- sdboot:

  Bootstrap standard error

- z0:

  Bias-correction constant

- a:

  Acceleration constant

- sdjack:

  Jackknife standard error (if available)

- B:

  Number of bootstrap replicates

- boot_mean:

  Mean of bootstrap replicates

## Examples

``` r
data(diabetes, package = "bcaboot")
Xy <- cbind(diabetes$x, diabetes$y)
rfun <- function(Xy) {
  y <- Xy[, 11]; X <- Xy[, 1:10]
  summary(lm(y ~ X))$adj.r.squared
}
set.seed(1234)
result <- bca_nonpar(Xy, 1000, rfun, n_groups = 34, verbose = FALSE)
glance(result)
#> # A tibble: 1 × 9
#>   method accel      theta sdboot     z0       a sdjack     B boot_mean
#>   <chr>  <chr>      <dbl>  <dbl>  <dbl>   <dbl>  <dbl> <dbl>     <dbl>
#> 1 nonpar regression 0.507 0.0299 -0.332 -0.0127 0.0305  1000     0.517
```
