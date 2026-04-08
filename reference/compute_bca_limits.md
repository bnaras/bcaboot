# Compute BCa confidence limits from bootstrap replicates.

BCa percentile formula: maps nominal alpha levels to
bias-and-acceleration adjusted percentiles of the bootstrap
distribution. Efron (1987) Sec 2.

## Usage

``` r
compute_bca_limits(z0, a, zalpha, boot_reps)
```

## Arguments

- z0:

  Bias-correction constant (scalar).

- a:

  Acceleration constant (scalar).

- zalpha:

  Normal quantiles at desired alpha levels (vector).

- boot_reps:

  Sorted vector of B bootstrap replicates.

## Value

Named list with `limits` (BCa quantiles), `bca_idx` (indices into sorted
boot_reps), and `iles` (adjusted percentile levels).
