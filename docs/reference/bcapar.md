# Compute parametric bootstrap confidence intervals

bcapar computes parametric bootstrap confidence intervals for a
real-valued parameter theta in a p-parameter exponential family. It is
described in Section 4 of the reference below.

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

  A vector of parametric bootstrap replications of theta of length `B`,
  usually large, say `B = 2000`

- bb:

  A `B` by `p` matrix of natural sufficient vectors, where `p` is the
  dimension of the exponential family.

- alpha:

  percentiles desired for the bca confidence limits. One only needs to
  provide `alpha` values below 0.5; the upper limits are automatically
  computed

- J, K:

  Parameters controlling the jackknife estimates of Monte Carlo error:
  `J` jackknife folds, with the jackknife standard errors averaged over
  `K` random divisions of `bb`

- trun:

  Truncation parameter used in the calculation of the acceleration `a`.

- pct:

  Proportion of "nearby" b vectors used in the calculation of `t.`, the
  gradient vector of theta.

- cd:

  If cd is 1 the bca confidence density is also returned; see Section
  11.6 in reference Efron and Hastie (2016) below

- func:

  Function \\\hat{\theta} = func(b)\\. If this is not missing then
  output includes *abc* estimates; see reference DiCiccio and
  Efron (1992) below

## Value

a named list of several items:

- **lims** : Bca confidence limits (first column) and the standard
  limits (fourth column). Also the abc limits (fifth column) if `func`
  is provided. The second column, `jacksd`, are the jackknife estimates
  of Monte Carlo error; `pct`, the third column are the proportion of
  the replicates `tt` less than each `bcalim` value

- **stats** : Estimates and their jackknife Monte Carlo errors: `theta`
  = \\\hat{\theta}\\; `sd`, the bootstrap standard deviation for
  \\\hat{\theta}\\; `a` the acceleration estimate; `az` another
  acceleration estimate that depends less on extreme values of `tt`;
  `z0` the bias-correction estimate; `A` the big-A measure of raw
  acceleration; `sdd` delta method estimate for standard deviation of
  \\\hat{\theta}\\; `mean` the average of `tt`

- **abcstats** : The abc estimates of `a` and `z0`, returned if `func`
  was provided

- **ustats** : The bias-corrected estimator `2 * t0 - mean(tt)`.
  `ustats` gives `ustat`, an estimate `sdu` of its sampling error, and
  jackknife estimates of monte carlo error for both `ustat` and `sdu`.
  Also given is `B`, the number of bootstrap replications

- **seed** : The random number state for reproducibility

## References

DiCiccio T and Efron B (1996). Bootstrap confidence intervals.
Statistical Science 11, 189-228

T. DiCiccio and B. Efron. More accurate confidence intervals in
exponential families. Biometrika (1992) p231-245.

Efron B (1987). Better bootstrap confidence intervals. JASA 82, 171-200

B. Efron and T. Hastie. Computer Age Statistical Inference. Cambridge
University Press, 2016.

B. Efron and B. Narasimhan. Automatic Construction of Bootstrap
Confidence Intervals, 2018.

## Examples

``` r
data(diabetes, package = "bcaboot")
X <- diabetes$x
y <- scale(diabetes$y, center = TRUE, scale = FALSE)
lm.model <- lm(y ~ X - 1)
mu.hat <- lm.model$fitted.values
sigma.hat <- stats::sd(lm.model$residuals)
t0 <- summary(lm.model)$adj.r.squared
y.star <- sapply(mu.hat, rnorm, n = 1000, sd = sigma.hat)
tt <- apply(y.star, 1, function(y) summary(lm(y ~ X - 1))$adj.r.squared)
b.star <- y.star %*% X
set.seed(1234)
bcapar(t0 = t0, tt = tt, bb = b.star)
#> Warning: `bcapar()` was deprecated in bcaboot 1.0.
#> ℹ Please use `bca_par()` instead.
#> 
#> ── BCa Bootstrap Confidence Intervals 
#> Method: par
#> B = 1000, theta = 0.5065862, sdboot = 0.02803423
#> 
#> Confidence limits:
#>  conf.level    bca.lo    bca.hi    std.lo    std.hi
#>        0.95 0.4381990 0.5514613 0.4516401 0.5615323
#>        0.90 0.4450255 0.5423186 0.4604740 0.5526984
#>        0.80 0.4563314 0.5317548 0.4706589 0.5425135
#>        0.68 0.4672079 0.5248057 0.4787074 0.5344651
#> 
#> Diagnostics:
#> z0 = -0.3584588, a = 0.0019184
```
