# Nonparametric bias-corrected and accelerated bootstrap confidence limits

This function is a version of `bcajack` that allows all the
recomputations of the original statistic function \\f\\ to be carried
out separately. This is an advantage if \\f\\ is time-consuming, in
which case the B replications for the nonparametric bca calculations
might need to be done on a distributed basis.

To use `bcajack2` in this mode, we first compute a list `Blist` via
`Blist <- list(Y = Y, tt = tt, t0 = t0)`. Here `tt` is a vector of
length `B` having i-th entry `tt[i] <- func(x[Ii,], ...)`, where `x` is
the \\n \times p\\ data matrix and `Ii` is a bootstrap vector of
(observation) indices. `Y` is a `B` by \\n\\ count matrix, whose i-th
row is the counts corresponding to `Ii`. For example if n = 5 and
`Ii = (2, 5, 2, 1, 4)`, then `Yi = (1, 2, 0, 1, 1)`. Having computed
`Blist`, `bcajack2` is invoked as `bcajack2(Blist)` without need to
enter the function \\func\\.

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

  an \\n \times p\\ data matrix, rows are observed \\p\\-vectors,
  assumed to be independently sampled from target population. If \\p\\
  is 1 then `x` can be a vector.

- B:

  number of bootstrap replications. `B` can also be a vector of `B`
  bootstrap replications of the estimated parameter of interest,
  computed separately. If `B` is `Blist` as explained above, `x` is not
  needed.

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

- pct:

  `bcajack2` uses those count vectors nearest (1,1,...1) to estimate the
  gradient of the statistic, "nearest" being defined as those count
  vectors in the smallest `pct` of all B of them. Default value for
  \`pct is 1/3 (see appendix in Efron and Narasimhan for further
  details)

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

- **stats** : top line of stats shows 5 estimates: theta is \\func(x)\\,
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
bcajack2(x = Xy, B = 1000, func = rfun, m = 40, verbose = FALSE)
#> Warning: `bcajack2()` was deprecated in bcaboot 1.0.
#> ℹ Please use `bca_nonpar()` instead.
#> BCa Bootstrap Confidence Intervals
#>   Method: nonpar (regression acceleration)
#>   B = 1000, theta = 0.5065603, sdboot = 0.03076788
#> 
#> Confidence limits:
#>  conf.level    bca.lo    bca.hi    std.lo    std.hi
#>        0.95 0.4289798 0.5540292 0.4462564 0.5668642
#>        0.90 0.4374506 0.5459683 0.4559517 0.5571690
#>        0.80 0.4533929 0.5357038 0.4671297 0.5459909
#>        0.68 0.4621266 0.5255210 0.4759630 0.5371577
#> 
#> Diagnostics:
#>   z0 = -0.3853205, a = -0.01094242, sdjack = 0.03139928
```
