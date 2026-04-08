# bcaboot 1.0

## Breaking changes

* New primary API: `bca_nonpar()` replaces `bcajack()`, `bcajack2()`,
  and `bcanon()`. `bca_par()` replaces `bcapar()`.

* All computation functions now return a canonical `"bcaboot"` object
  with consistent structure: `$limits` (9x4 matrix), `$stats` (2x5
  matrix), `$B_mean`, `$ustats`. The old field names (`$lims`,
  `$B.mean`) are no longer used.

* User-facing parameter `conf.level` (e.g., `c(0.95, 0.90)`) replaces
  `alpha` (e.g., `c(0.025, 0.05)`). Internally, alpha expansion is
  unchanged.

* Descriptive parameter names: `n_groups` (was `m`), `group_reps`
  (was `mr`), `kl_fraction` (was `pct`), `n_jack` (was `K`),
  `jack_groups` (was `J`), `truncation` (was `trun`),
  `conf_density` (was `cd`), `boot_data` (was overloaded `B`).

* `bcanon()` removed (was never on CRAN).

## Deprecations

* `bcajack()`, `bcajack2()`, `bcapar()`, and `bcaplot()` are
 deprecated with `lifecycle` warnings (once per session). They will
  continue to work through thin wrappers that translate to the new
  API.

## New features
 
* `tidy()` method: returns a tibble with one row per (confidence
  level, method) pair, following broom conventions (`conf.level`,
  `method`, `estimate`, `conf.low`, `conf.high`).

* `glance()` method: returns a one-row tibble summarizing the
  bootstrap run (`method`, `accel`, `theta`, `sdboot`, `z0`, `a`,
  `sdjack`, `B`, `boot_mean`).

* `autoplot()` method: ggplot2 visualization of BCa vs standard
  confidence intervals. Requires ggplot2 (in Suggests).

* Improved `print()`: formatted summary with confidence limits table
  instead of raw list dump.

* `bca_nonpar()` always computes the gbca diagnostic when using
  `accel = "regression"` (with `tryCatch` for robustness).

## Dependencies

* Added `tibble`, `generics`, `lifecycle` to Imports.
* Added `ggplot2` to Suggests (for `autoplot()`).

## Internal

* Shared helpers extracted to `R/bca_core.R`: `expand_alpha()`,
  `bootstrap_resample()`, `jackknife_accel()`, `regression_accel()`,
  `compute_ustats()`.
* `new_bcaboot()` canonical constructor guarantees consistent return
  structure across all computation functions.
* Test suite expanded from 14 to 34 tests, including adapter tests
  proving numerical equivalence against pre-1.0 fixtures.
