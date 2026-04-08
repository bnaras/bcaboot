# Nonparametric BCa bootstrap confidence intervals

Computes nonparametric bias-corrected and accelerated (BCa) bootstrap
confidence intervals. Consolidates the functionality of the deprecated
`bcajack`, `bcajack2`, and `bcanon` functions.

## Usage

``` r
bca_nonpar(
  x,
  B,
  func,
  ...,
  accel = c("regression", "jackknife"),
  conf.level = c(0.95, 0.9, 0.8, 0.68),
  n_groups = nrow(x),
  group_reps = 5,
  kl_fraction = 0.333,
  n_jack = 2,
  jack_groups = 12,
  boot_data = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  An \\n \times p\\ data matrix. Rows are observed \\p\\-vectors,
  assumed independently sampled. If \\p = 1\\ then `x` can be a vector.

- B:

  Number of bootstrap replications (integer).

- func:

  Function computing the parameter of interest. `func(x, ...)` should
  return a scalar for any data matrix.

- ...:

  Additional arguments passed to `func`.

- accel:

  Method for estimating acceleration: `"regression"` (default) uses
  local regression on bootstrap count vectors nearest to uniform;
  `"jackknife"` uses classical delete-one (or delete-group) jackknife
  influence values.

- conf.level:

  Confidence levels for the intervals (e.g.,
  `c(0.95, 0.90, 0.80, 0.68)`). Only values in `(0, 1)` are used.

- n_groups:

  Number of groups for jackknife. Default `nrow(x)` gives delete-one
  jackknife. Set smaller (e.g., 20-40) to speed up computation on large
  datasets.

- group_reps:

  When `n_groups < nrow(x)`, number of random grouping repetitions to
  average over. Only used with `accel = "jackknife"`.

- kl_fraction:

  Fraction of bootstrap count vectors nearest to uniform used in the
  local regression. Only used with `accel = "regression"`. Default 1/3.

- n_jack:

  Number of delete-d jackknife repetitions for estimating internal
  (Monte Carlo) standard error. Set to 0 to skip.

- jack_groups:

  Number of groups per jackknife fold.

- boot_data:

  Optional pre-computed bootstrap data: a list with components `Y` (B x
  n count matrix), `tt` (bootstrap replicates), `t0` (original
  estimate). When provided, `B` is ignored and no resampling is
  performed.

- verbose:

  Logical; show progress bar during bootstrap.

## Value

An object of class `"bcaboot"` with components:

- limits:

  9-row matrix of confidence limits (columns: `bca`, `jacksd`, `std`,
  `pct`)

- stats:

  2-row matrix of estimates and jackknife SEs (columns: `theta`,
  `sdboot`, `z0`, `a`, `sdjack`)

- B_mean:

  Bootstrap sample size and mean of replicates

- ustats:

  Bias-corrected estimator and its SE

- diagnostic:

  gbca diagnostic (when `accel = "regression"`)

## References

Efron B (1987). Better bootstrap confidence intervals. JASA 82, 171-200.

Efron B and Narasimhan B (2020). Automatic construction of bootstrap
confidence intervals. JASA 115, 1895-1905.

## Examples

``` r
data(diabetes, package = "bcaboot")
Xy <- cbind(diabetes$x, diabetes$y)
rfun <- function(Xy) {
  y <- Xy[, 11]
  X <- Xy[, 1:10]
  summary(lm(y ~ X))$adj.r.squared
}
set.seed(1234)
result <- bca_nonpar(x = Xy, B = 1000, func = rfun, n_groups = 34,
                     verbose = FALSE)
print(result)
#> 
#> ── BCa Bootstrap Confidence Intervals 
#> Method: nonpar (regression acceleration)
#> B = 1000, theta = 0.5065603, sdboot = 0.02988642
#> 
#> Confidence limits:
#>  conf.level    bca.lo    bca.hi    std.lo    std.hi
#>        0.95 0.4205204 0.5544006 0.4479840 0.5651366
#>        0.90 0.4473162 0.5458156 0.4574015 0.5557191
#>        0.80 0.4574096 0.5349201 0.4682593 0.5448613
#>        0.68 0.4654674 0.5272497 0.4768395 0.5362811
#> 
#> Diagnostics:
#> z0 = -0.3318533, a = -0.01267648, sdjack = 0.03047356
```
