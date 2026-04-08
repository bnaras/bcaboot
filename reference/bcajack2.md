# Nonparametric BCa via regression (deprecated)

**\[deprecated\]**

`bcajack2()` was deprecated in bcaboot 1.0 in favour of
[`bca_nonpar()`](bca_nonpar.md).

## Usage

``` r
bcajack2(
  x,
  B,
  func,
  ...,
  m = NULL,
  mr,
  pct = 0.333,
  K = 2,
  J = 12,
  alpha = c(0.025, 0.05, 0.1, 0.16),
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

- m:

  Number of groups (now `n_groups`).

- mr:

  Unused (included for compatibility).

- pct:

  KL fraction (now `kl_fraction`).

- K:

  Jackknife repetitions (now `n_jack`).

- J:

  Jackknife groups (now `jack_groups`).

- alpha:

  Alpha levels (now use `conf.level`).

- verbose:

  Logical; show progress bar during bootstrap.
