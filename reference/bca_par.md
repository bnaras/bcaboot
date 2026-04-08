# Parametric BCa bootstrap confidence intervals

Computes parametric bootstrap confidence intervals for a real-valued
parameter theta in a p-parameter exponential family.

## Usage

``` r
bca_par(
  t0,
  tt,
  bb,
  conf.level = c(0.95, 0.9, 0.8, 0.68),
  n_jack = 6,
  jack_groups = 10,
  truncation = 0.001,
  kl_fraction = 0.333,
  conf_density = FALSE,
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

- conf.level:

  Confidence levels for the intervals.

- n_jack:

  Number of jackknife repetitions for internal SE. Set to 0 to skip.

- jack_groups:

  Number of groups per jackknife fold.

- truncation:

  Truncation parameter for acceleration calculation. Can be a vector to
  compute acceleration at multiple levels.

- kl_fraction:

  Proportion of "nearby" b vectors used in gradient estimation.

- conf_density:

  Logical; if TRUE, return BCa confidence density weights.

- func:

  Optional function for ABC (analytical bootstrap) limits.

## Value

An object of class `"bcaboot"` with components:

- limits:

  9-row matrix of confidence limits

- stats:

  2-row matrix of estimates and jackknife SEs

- B_mean:

  Bootstrap sample size and mean of replicates

- ustats:

  Bias-corrected estimator and its SE

- abc:

  ABC limits and stats (when `func` provided)

- conf_density:

  Confidence density weights (when requested)

- accel_matrix:

  Acceleration at multiple truncation levels

## References

Efron B (1987). Better bootstrap confidence intervals. JASA 82, 171-200.

DiCiccio T and Efron B (1992). More accurate confidence intervals in
exponential families. Biometrika, 231-245.

Efron B and Narasimhan B (2020). The Automatic Construction of Bootstrap
Confidence Intervals. Journal of Computational and Graphical Statistics,
29(3), 608-619.
[doi:10.1080/10618600.2020.1714633](https://doi.org/10.1080/10618600.2020.1714633)

## Examples

``` r
data(diabetes, package = "bcaboot")
X <- diabetes$x
y <- scale(diabetes$y, center = TRUE, scale = FALSE)
lm.model <- lm(y ~ X - 1)
mu.hat <- lm.model$fitted.values
sigma.hat <- stats::sd(lm.model$residuals)
t0 <- summary(lm.model)$adj.r.squared
y.star <- sapply(mu.hat, rnorm, n = 500, sd = sigma.hat)
tt <- apply(y.star, 1, function(y) summary(lm(y ~ X - 1))$adj.r.squared)
b.star <- y.star %*% X
set.seed(1234)
bca_par(t0 = t0, tt = tt, bb = b.star)
#> 
#> ── BCa Bootstrap Confidence Intervals 
#> Method: par
#> B = 500, theta = 0.5065862, sdboot = 0.0297771
#> 
#> Confidence limits:
#>  conf.level    bca.lo    bca.hi    std.lo    std.hi
#>        0.95 0.4472066 0.5565305 0.4482242 0.5649483
#>        0.90 0.4518055 0.5478429 0.4576072 0.5555652
#>        0.80 0.4626100 0.5349231 0.4684253 0.5447471
#>        0.68 0.4706893 0.5267426 0.4769741 0.5361983
#> 
#> Diagnostics:
#> z0 = -0.3107377, a = 0.01129753
```
