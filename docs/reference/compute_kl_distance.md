# Compute Poisson/multinomial deviance from uniform for bootstrap count matrix.

Each row of Y is a bootstrap count vector. The deviance measures how far
each bootstrap sample departs from simple random sampling (uniform
counts). Limit when Yi=0: the term equals 2. EN20 appendix.

## Usage

``` r
compute_kl_distance(Y)
```

## Arguments

- Y:

  B x m count matrix (rows = bootstrap samples, cols = observation
  counts).

## Value

Numeric vector of length B with deviance values.
