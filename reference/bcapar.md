# Parametric BCa (deprecated)

**\[deprecated\]**

`bcapar()` was deprecated in bcaboot 1.0 in favour of
[`bca_par()`](bca_par.md).

## Usage

``` r
bcapar(
  t0,
  tt,
  bb,
  alpha = c(0.025, 0.05, 0.1, 0.16),
  J = 10,
  K = 6,
  trun = 0.001,
  pct = 0.333,
  cd = 0,
  func
)
```

## Arguments

- t0:

  Observed estimate of theta, usually by maximum likelihood.

- tt:

  A vector of B parametric bootstrap replications of theta.

- bb:

  A B by p matrix of natural sufficient vectors, where p is the
  dimension of the exponential family.

- alpha:

  Alpha levels (now use `conf.level`).

- J:

  Jackknife groups (now `jack_groups`).

- K:

  Jackknife repetitions (now `n_jack`).

- trun:

  Truncation (now `truncation`).

- pct:

  KL fraction (now `kl_fraction`).

- cd:

  Confidence density flag (now `conf_density`).

- func:

  Optional function for ABC (analytical bootstrap) limits.
