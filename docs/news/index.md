# Changelog

## bcaboot 1.0

### Breaking changes

- New primary API: `bca_nonpar()` replaces
  [`bcajack()`](../reference/bcajack.md),
  [`bcajack2()`](../reference/bcajack2.md), and
  [`bcanon()`](../reference/bcanon.md). `bca_par()` replaces
  [`bcapar()`](../reference/bcapar.md).

- All computation functions now return a canonical `"bcaboot"` object
  with consistent structure: `$limits` (9x4 matrix), `$stats` (2x5
  matrix), `$B_mean`, `$ustats`. The old field names (`$lims`,
  `$B.mean`) are no longer used.

- User-facing parameter `conf.level` (e.g., `c(0.95, 0.90)`) replaces
  `alpha` (e.g., `c(0.025, 0.05)`). Internally, alpha expansion is
  unchanged.

- Descriptive parameter names: `n_groups` (was `m`), `group_reps` (was
  `mr`), `kl_fraction` (was `pct`), `n_jack` (was `K`), `jack_groups`
  (was `J`), `truncation` (was `trun`), `conf_density` (was `cd`),
  `boot_data` (was overloaded `B`).

- [`bcanon()`](../reference/bcanon.md) removed (was never on CRAN).

### Deprecations

- [`bcajack()`](../reference/bcajack.md),
  [`bcajack2()`](../reference/bcajack2.md),
  [`bcapar()`](../reference/bcapar.md), and
  [`bcaplot()`](../reference/bcaplot.md) are deprecated with `lifecycle`
  warnings (once per session). They will continue to work through thin
  wrappers that translate to the new API.

### Bug fixes

- **Issue [\#2](https://github.com/bnaras/bcaboot/issues/2)**: Fixed
  grouped jackknife (`n_groups < n`) using
  `sapply(seq_len_m, sample.int, ...)` which accidentally passed
  iteration values as the first positional arg, overriding the `n=`
  keyword and enabling replacement. Replaced with
  `matrix(sample.int(n, n-r), nrow=m)` for correct without-replacement
  partition. (Reported by admash)

- **Issue [\#7](https://github.com/bnaras/bcaboot/issues/7)**: Fixed
  dimension drop when `x` is a single-column matrix (vector input).
  `x[-i, ]` on an `n x 1` matrix returned a vector instead of a matrix,
  causing type inconsistency between jackknife and bootstrap code paths.
  Added `drop = FALSE` to all matrix subsetting passed to `func`.
  (Reported by R180)

- **Issue [\#4](https://github.com/bnaras/bcaboot/issues/4) /
  [\#1](https://github.com/bnaras/bcaboot/issues/1)**:
  `regression_accel` now checks for underdetermined regression
  (`nearby samples < ncol(Y)`) and errors with an actionable message
  instead of silently returning NAs. Suggests increasing `B`, increasing
  `kl_fraction`, or switching to `accel = "jackknife"`. (Reported by Tim
  Pollington, Thomas Covert)

- **PR [\#8](https://github.com/bnaras/bcaboot/issues/8)**: Fixed
  missing [`return()`](https://rdrr.io/r/base/function.html) in `K = 0`
  early-exit path of [`bcajack()`](../reference/bcajack.md),
  [`bcajack2()`](../reference/bcajack2.md), and
  [`bcapar()`](../reference/bcapar.md), which caused execution to fall
  through to the `K > 0` code. In the rewrite, all early exits use
  explicit `return(new_bcaboot(...))`. (Reported by Bettina Gruen)

### New features

- `tidy()` method: returns a tibble with one row per (confidence level,
  method) pair, following broom conventions (`conf.level`, `method`,
  `estimate`, `conf.low`, `conf.high`).

- `glance()` method: returns a one-row tibble summarizing the bootstrap
  run (`method`, `accel`, `theta`, `sdboot`, `z0`, `a`, `sdjack`, `B`,
  `boot_mean`).

- [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  method: ggplot2 visualization of BCa vs standard confidence intervals.
  Requires ggplot2 (in Suggests).

- Improved [`print()`](https://rdrr.io/r/base/print.html): formatted
  summary with confidence limits table instead of raw list dump.

- `bca_nonpar()` always computes the gbca diagnostic when using
  `accel = "regression"` (with `tryCatch` for robustness).

- All user-facing messages use `cli` for consistent formatting.

### Dependencies

- Added `cli`, `tibble`, `generics`, `lifecycle` to Imports.
- Added `ggplot2` to Suggests (for
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)).

### Internal

- Shared helpers extracted to `R/bca_core.R`: `expand_alpha()`,
  `bootstrap_resample()`, `jackknife_accel()`, `regression_accel()`,
  `compute_ustats()`.
- `new_bcaboot()` canonical constructor guarantees consistent return
  structure across all computation functions.
- Test suite expanded from 14 to 34 tests, including adapter tests
  proving numerical equivalence against pre-1.0 fixtures.
