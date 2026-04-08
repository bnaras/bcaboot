## Regression tests for bcaboot
## Adapter tests: call new API, remap to old structure, compare against
## saved baseline fixtures to prove numerical equivalence.

library(bcaboot)

fixture_dir <- system.file("tinytest", "fixtures", package = "bcaboot")
tol <- 1e-10

## ---------------------------------------------------------------
## Common setup: diabetes data and rfun
## ---------------------------------------------------------------
data(diabetes, package = "bcaboot")
Xy <- cbind(diabetes$x, diabetes$y)
rfun <- function(Xy) {
    y <- Xy[, 11]
    X <- Xy[, 1:10]
    summary(lm(y ~ X))$adj.r.squared
}

## ---------------------------------------------------------------
## bca_nonpar(accel="jackknife") vs old bcajack fixtures
## ---------------------------------------------------------------
ref <- readRDS(file.path(fixture_dir, "ref_bcajack.rds"))
set.seed(1234)
cur <- bca_nonpar(x = Xy, B = 1000, func = rfun,
                  accel = "jackknife", n_groups = 34,
                  n_jack = 2, jack_groups = 10,
                  verbose = FALSE)

## BCa and std limits must match exactly
expect_equal(cur$limits[, "bca"], ref$limits[, "bca"], tolerance = tol,
             info = "bcajack equiv: bca limits match")
expect_equal(cur$limits[, "std"], ref$limits[, "std"], tolerance = tol,
             info = "bcajack equiv: std limits match")
## Stats estimate row must match
expect_equal(unname(cur$stats["est", ]), unname(ref$stats["est", ]),
             tolerance = tol, info = "bcajack equiv: stats est match")
## Ustats must match
expect_equal(unname(cur$ustats), unname(ref$ustats), tolerance = tol,
             info = "bcajack equiv: ustats match")

## ---------------------------------------------------------------
## bca_nonpar(accel="regression") vs old bcajack2 fixtures
## ---------------------------------------------------------------
ref2 <- readRDS(file.path(fixture_dir, "ref_bcajack2.rds"))
set.seed(1234)
cur2 <- bca_nonpar(x = Xy, B = 1000, func = rfun,
                   accel = "regression", n_groups = 40,
                   n_jack = 2, jack_groups = 12,
                   verbose = FALSE)

expect_equal(cur2$limits[, "bca"], ref2$lims[, "bca"], tolerance = tol,
             info = "bcajack2 equiv: bca limits match")
expect_equal(cur2$limits[, "std"], ref2$lims[, "std"], tolerance = tol,
             info = "bcajack2 equiv: std limits match")
expect_equal(unname(cur2$stats["est", ]), unname(ref2$stats["est", ]),
             tolerance = tol, info = "bcajack2 equiv: stats est match")
expect_equal(unname(cur2$ustats), unname(ref2$ustats), tolerance = tol,
             info = "bcajack2 equiv: ustats match")

## ---------------------------------------------------------------
## bca_nonpar(accel="regression") vs old bcanon fixtures
## (bcanon uses regression + gbca diagnostic)
## ---------------------------------------------------------------
d <- readRDS(file.path(fixture_dir, "bcanon_data.rds"))
func5 <- function(Xy) {
    y <- Xy[, 7]; X <- Xy[, 1:6]
    re <- glm(y ~ X, binomial)
    re$coef[4]
}
refn <- readRDS(file.path(fixture_dir, "ref_bcanon.rds"))
set.seed(1234)
curn <- bca_nonpar(x = d$Xy, B = 1000, func = func5,
                   accel = "regression",
                   n_jack = 2, jack_groups = 12,
                   verbose = FALSE)

expect_equal(curn$limits[, "bca"], refn$lims[, "bca"], tolerance = tol,
             info = "bcanon equiv: bca limits match")
expect_equal(curn$limits[, "std"], refn$lims[, "std"], tolerance = tol,
             info = "bcanon equiv: std limits match")
## bcanon stats has different columns (theta,sdboot,sdjack,z0,a,A,sese)
## but the common 5 must match
expect_equal(unname(curn$stats["est", "theta"]),
             unname(refn$stats["estimate", "theta"]), tolerance = tol,
             info = "bcanon equiv: theta match")
expect_equal(unname(curn$stats["est", "z0"]),
             unname(refn$stats["estimate", "z0"]), tolerance = tol,
             info = "bcanon equiv: z0 match")
expect_equal(unname(curn$stats["est", "a"]),
             unname(refn$stats["estimate", "a"]), tolerance = tol,
             info = "bcanon equiv: a match")

## gbca diagnostic present
expect_true(!is.null(curn$diagnostic),
            info = "bcanon equiv: diagnostic present")

## ---------------------------------------------------------------
## bca_par vs old bcapar fixtures
## ---------------------------------------------------------------
X <- diabetes$x
y <- scale(diabetes$y, center = TRUE, scale = FALSE)
lm_model <- lm(y ~ X - 1)
mu_hat <- lm_model$fitted.values
sigma_hat <- stats::sd(lm_model$residuals)
t0 <- summary(lm_model)$adj.r.squared
set.seed(1234)
y_star <- sapply(mu_hat, rnorm, n = 1000, sd = sigma_hat)
tt <- apply(y_star, 1, function(y) summary(lm(y ~ X - 1))$adj.r.squared)
b_star <- y_star %*% X

refp <- readRDS(file.path(fixture_dir, "ref_bcapar.rds"))
set.seed(5678)
curp <- bca_par(t0 = t0, tt = tt, bb = b_star)

expect_equal(curp$limits[, "bca"], refp$lims[, "bca"], tolerance = tol,
             info = "bcapar equiv: bca limits match")
expect_equal(curp$limits[, "std"], refp$lims[, "std"], tolerance = tol,
             info = "bcapar equiv: std limits match")

## ---------------------------------------------------------------
## Structure tests: canonical return format
## ---------------------------------------------------------------
set.seed(42)
r <- bca_nonpar(Xy, 500, rfun, n_groups = 34, n_jack = 0, verbose = FALSE)

expect_true(inherits(r, "bcaboot"), info = "class is bcaboot")
expect_equal(dim(r$limits), c(9L, 4L), info = "limits is 9x4")
expect_equal(colnames(r$limits), c("bca", "jacksd", "std", "pct"),
             info = "limits colnames")
expect_equal(dim(r$stats), c(2L, 5L), info = "stats is 2x5")
expect_equal(rownames(r$stats), c("est", "jsd"), info = "stats rownames")
expect_equal(colnames(r$stats), c("theta", "sdboot", "z0", "a", "sdjack"),
             info = "stats colnames")
expect_equal(length(r$B_mean), 2L, info = "B_mean length 2")
expect_equal(length(r$ustats), 2L, info = "ustats length 2")
expect_equal(r$method, "nonpar", info = "method is nonpar")

## ---------------------------------------------------------------
## tidy() output tests
## ---------------------------------------------------------------
td <- tidy(r)
expect_true(inherits(td, "tbl_df"), info = "tidy returns tibble")
expect_equal(ncol(td), 7L, info = "tidy has 7 columns")
expect_true(all(c("conf.level", "method", "estimate", "conf.low",
                  "conf.high") %in% names(td)),
            info = "tidy has required columns")
expect_equal(nrow(td), 8L, info = "tidy has 8 rows (4 coverage x 2 methods)")
expect_true(all(td$method %in% c("bca", "standard")),
            info = "tidy methods are bca and standard")

## ---------------------------------------------------------------
## glance() output tests
## ---------------------------------------------------------------
gl <- glance(r)
expect_true(inherits(gl, "tbl_df"), info = "glance returns tibble")
expect_equal(nrow(gl), 1L, info = "glance has 1 row")
expect_true(all(c("method", "accel", "theta", "sdboot", "z0", "a",
                  "B") %in% names(gl)),
            info = "glance has required columns")

## ---------------------------------------------------------------
## autoplot test (skip if ggplot2 not available)
## ---------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::autoplot(r)
    expect_true(inherits(p, "gg"), info = "autoplot returns ggplot object")
}
