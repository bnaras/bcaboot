# Compute nonparametric BCa bootstrap confidence intervals

This function computes nonparametric BCa confidence intervals using
regression on bootstrap count vectors (same approach as
[bcajack2](bcajack2.md)) plus a generalized BCa (gbca) diagnostic that
tests adequacy of the BCa approximation via warped normal comparison.
Efron (1987) Sec 7.

## Usage

``` r
bcanon(
  B,
  x,
  func,
  ...,
  m = nrow(x),
  pct = 0.333,
  K = 2,
  J = 12,
  alpha = c(0.025, 0.05, 0.1, 0.16),
  verbose = TRUE
)
```

## Arguments

- B:

  number of bootstrap replications, or a list with components `Y` (count
  matrix), `tt` (bootstrap replicates), `t0` (original estimate)

- x:

  data matrix of dimension \\n\times p\\, rows assumed mutually
  independent

- func:

  R function, where \\func(x, \ldots)\\ is the statistic of interest

- ...:

  additional arguments for `func`

- m:

  number of units; if \\m\<n\\ then original units collected into groups
  each of size \\n/m\\; useful if \\n\\ is very large

- pct:

  proportion of nearby count vectors used in finding gradient

- K:

  number of jackknife repetitions for internal standard error estimation

- J:

  number of jackknife folds per repetition

- alpha:

  percentiles of coverage probabilities for bootstrap intervals

- verbose:

  logical for verbose progress messages

## Value

a named list of class `bcaboot` containing:

- **call**: the matched call

- **lims** : BCa confidence limits (`bca`), jackknife internal standard
  errors (`jacksd`), standard limits (`std`), and bootstrap percentiles
  (`pct`)

- **stats** : Estimates and jackknife Monte Carlo errors for `theta`,
  `sdboot`, `sdjack`, `z0`, `a`, plus big-A acceleration and SE
  stability (`A`, `sese`)

- **B.mean** : bootstrap sample size B and mean of replications

- **equiv** : gbca coverage estimates showing `alpha` vs `alphatilda`
  (coverage-adjusted alpha); if BCa is exact these are equal. Also
  jackknife Monte Carlo standard deviations.

- **seed** : the random number state for reproducibility
