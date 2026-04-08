# Nonparametric bias-corrected and accelerated bootstrap confidence limits

This routine computes nonparametric confidence intervals for bootstrap
estimates. For reproducibility, save or set the random number state
before calling this routine.

## Usage

``` r
bcajack(
  x,
  B,
  func,
  ...,
  m = nrow(x),
  mr = 5,
  K = 2,
  J = 10,
  alpha = c(0.025, 0.05, 0.1, 0.16),
  verbose = TRUE
)
```

## Arguments

- x:

  an \\n \times p\\ data matrix, rows are observed \\p\\-vectors,
  assumed to be independently sampled from target population. If \\p\\
  is 1 then `x` can be a vector.

- B:

  number of bootstrap replications. It can also be a vector of `B`
  bootstrap replications of the estimated parameter of interest,
  computed separately.

- func:

  function \\\hat{\theta}=func(x)\\ computing estimate of the parameter
  of interest; \\func(x)\\ should return a real value for any \\n^\prime
  \times p\\ matrix \\x^\prime\\, \\n^\prime\\ not necessarily equal to
  \\n\\

- ...:

  additional arguments for `func`.

- m:

  an integer less than or equal to \\n\\; the routine collects the \\n\\
  rows of `x` into `m` groups to speed up the jackknife calculations for
  estimating the acceleration value \\a\\; typically `m` is 20 or 40 and
  does not have to exactly divide \\n\\. However, warnings will be
  shown.

- mr:

  if \\m \< n\\ then `mr` repetions of the randomly grouped jackknife
  calculations are averaged.

- K:

  a non-negative integer. If `K` \> 0, bcajack also returns estimates of
  *internal standard error*, that is, of the variability due to stopping
  at `B` bootstrap replications rather than going on to infinity. These
  are obtained from a second type of jackknifing, taking an average of
  `K` separate jackknife estimates, each randomly splitting the `B`
  bootstrap replications into `J` groups.

- J:

  the number of groups into which the bootstrap replications are split

- alpha:

  percentiles desired for the bca confidence limits. One only needs to
  provide `alpha` values below 0.5; the upper limits are automatically
  computed

- verbose:

  logical for verbose progress messages

## Value

a named list of several items

- **lims** : first column shows the estimated bca confidence limits at
  the requested alpha percentiles. These can be compared with the
  standard limits \\\hat{\theta} + \hat{\sigma}z\_{\alpha}\\, third
  column. The second column `jacksd` gives the internal standard errors
  for the bca limits, quite small in the example. Column 4, `pct`, gives
  the percentiles of the ordered B bootstrap replications corresponding
  to the bca limits, eg the 897th largest replication equalling the .975
  bca limit .557.

- **stats** : top line of stats shows 5 estimates: theta is \\f(x)\\,
  original point estimate of the parameter of interest; `sdboot` is its
  bootstrap estimate of standard error; `z0` is the bca bias correction
  value, in this case quite negative; `a` is the *acceleration*, a
  component of the bca limits (nearly zero here); `sdjack` is the
  jackknife estimate of standard error for theta. Bottom line gives the
  internal standard errors for the five quantities above. This is
  substantial for `z0` above.

- **B.mean** : bootstrap sample size B, and the mean of the B bootstrap
  replications \\\hat{\theta^\*}\\

- **ustats** : The bias-corrected estimator `2 * t0 - mean(tt)`, and an
  estimate `sdu` of its sampling error

- **seed** : The random number state for reproducibility

## Details

Bootstrap confidence intervals depend on three elements:

- the cdf of the \\B\\ bootstrap replications \\t_i^\*\\, \\i=1\ldots
  B\\

- the bias-correction number \\z_0=\Phi(\sum_i^B I(t_i^\* \< t_0) / B
  )\\ where \\t_0=f(x)\\ is the original estimate

- the acceleration number \\a\\ that measures the rate of change in
  \\\sigma\_{t_0}\\ as \\x\\, the data changes.

The first two of these depend only on the bootstrap distribution, and
not how it is generated: parametrically or non-parametrically. Program
bcajack can be used in a hybrid fashion in which the vector `tt` of B
bootstrap replications is first generated from a parametric model.

So, in the diabetes example below, we might first draw bootstrap samples
\\y^\* \sim N(X\hat{\beta}, \hat{\sigma}^2 I)\\ where \\\hat{\beta}\\
and \\\hat{\sigma}\\ were obtained from `lm(y~X)`; each \\y^\*\\ would
then provide a bootstrap replication `tstar = rfun(cbind(X, ystar))`.
Then we could get bca intervals from `bcajack(Xy, tt, rfun ....)` with
`tt`, the vector of B `tstar` values. The only difference from a full
parametric bca analysis would lie in the nonparametric estimation of
\\a\\, often a negligible error.

## References

DiCiccio T and Efron B (1996). Bootstrap confidence intervals.
Statistical Science 11, 189-228

Efron B (1987). Better bootstrap confidence intervals. JASA 82 171-200

B. Efron and B. Narasimhan. Automatic Construction of Bootstrap
Confidence Intervals, 2018.

## Examples

``` r
data(diabetes, package = "bcaboot")
Xy <- cbind(diabetes$x, diabetes$y)
rfun <- function(Xy) {
  y <- Xy[, 11]
  X <- Xy[, 1:10]
  summary(lm(y~X) )$adj.r.squared
}
set.seed(1234)
## n = 442 = 34 * 13
bcajack(x = Xy, B = 1000, func = rfun, m = 34, verbose = FALSE)
#> Warning: `bcajack()` was deprecated in bcaboot 1.0.
#> ℹ Please use `bca_nonpar()` instead.
#> 
#> ── BCa Bootstrap Confidence Intervals 
#> Method: nonpar (jackknife acceleration)
#> B = 1000, theta = 0.5065603, sdboot = 0.03284722
#> 
#> Confidence limits:
#>  conf.level    bca.lo    bca.hi    std.lo    std.hi
#>        0.95 0.4238753 0.5582875 0.4421809 0.5709397
#>        0.90 0.4342607 0.5496734 0.4525314 0.5605892
#>        0.80 0.4513913 0.5388794 0.4644649 0.5486557
#>        0.68 0.4617252 0.5305165 0.4738951 0.5392255
#> 
#> Diagnostics:
#> z0 = -0.26112, a = -0.01143494, sdjack = 0.02985439
```
