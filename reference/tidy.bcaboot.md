# Tidy a bcaboot result

Returns a tibble with one row per (confidence level, method)
combination, following broom conventions.

## Usage

``` r
# S3 method for class 'bcaboot'
tidy(x, ...)
```

## Arguments

- x:

  A `bcaboot` object from [`bca_nonpar()`](bca_nonpar.md) or
  [`bca_par()`](bca_par.md).

- ...:

  Ignored.

## Value

A tibble with columns:

- conf.level:

  Confidence level (e.g. 0.95)

- method:

  `"bca"` or `"standard"` (or `"abc"` when available)

- estimate:

  Point estimate (theta)

- conf.low:

  Lower confidence limit

- conf.high:

  Upper confidence limit

- jacksd.low:

  Jackknife SE of lower limit (BCa only; `NA` otherwise)

- jacksd.high:

  Jackknife SE of upper limit (BCa only; `NA` otherwise)

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
tidy(result)
#> # A tibble: 8 × 7
#>   conf.level method   estimate conf.low conf.high jacksd.low jacksd.high
#>        <dbl> <chr>       <dbl>    <dbl>     <dbl>      <dbl>       <dbl>
#> 1       0.95 bca         0.507    0.421     0.554    0.0106      0.00304
#> 2       0.95 standard    0.507    0.448     0.565   NA          NA      
#> 3       0.9  bca         0.507    0.447     0.546    0.00187     0.00176
#> 4       0.9  standard    0.507    0.457     0.556   NA          NA      
#> 5       0.8  bca         0.507    0.457     0.535    0.00284     0.00195
#> 6       0.8  standard    0.507    0.468     0.545   NA          NA      
#> 7       0.68 bca         0.507    0.465     0.527    0.00276     0.00180
#> 8       0.68 standard    0.507    0.477     0.536   NA          NA      
```
