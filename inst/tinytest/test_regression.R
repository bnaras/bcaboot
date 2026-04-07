## Regression tests for bcaboot
## Compare current output against saved baseline fixtures.

library(bcaboot)

fixture_dir <- system.file("tinytest", "fixtures", package = "bcaboot")
tol <- 1e-10

## ---------------------------------------------------------------
## Helper: strip seed from result (large, not meaningful to compare)
## ---------------------------------------------------------------
strip_seed <- function(x) { x$seed <- NULL; x }

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
## bcajack
## ---------------------------------------------------------------
ref <- readRDS(file.path(fixture_dir, "ref_bcajack.rds"))
set.seed(1234)
cur <- bcajack(x = Xy, B = 1000, func = rfun, m = 34, verbose = FALSE)

expect_equal(cur$lims, ref$lims, tolerance = tol,
             info = "bcajack: lims match baseline")
expect_equal(cur$stats, ref$stats, tolerance = tol,
             info = "bcajack: stats match baseline")
expect_equal(cur$B.mean, ref$B.mean, tolerance = tol,
             info = "bcajack: B.mean match baseline")
expect_equal(cur$ustats, ref$ustats, tolerance = tol,
             info = "bcajack: ustats match baseline")

## ---------------------------------------------------------------
## bcajack2
## ---------------------------------------------------------------
ref2 <- readRDS(file.path(fixture_dir, "ref_bcajack2.rds"))
set.seed(1234)
cur2 <- bcajack2(x = Xy, B = 1000, func = rfun, m = 40, verbose = FALSE)

expect_equal(cur2$lims, ref2$lims, tolerance = tol,
             info = "bcajack2: lims match baseline")
expect_equal(cur2$stats, ref2$stats, tolerance = tol,
             info = "bcajack2: stats match baseline")
expect_equal(cur2$B.mean, ref2$B.mean, tolerance = tol,
             info = "bcajack2: B.mean match baseline")
expect_equal(cur2$ustats, ref2$ustats, tolerance = tol,
             info = "bcajack2: ustats match baseline")

## ---------------------------------------------------------------
## bcanon
## ---------------------------------------------------------------
d <- readRDS(file.path(fixture_dir, "bcanon_data.rds"))
func5 <- function(Xy) {
    y <- Xy[, 7]; X <- Xy[, 1:6]
    re <- glm(y ~ X, binomial)
    re$coef[4]
}
refn <- readRDS(file.path(fixture_dir, "ref_bcanon.rds"))
set.seed(1234)
curn <- bcanon(B = 1000, x = d$Xy, func = func5, verbose = FALSE)

expect_equal(curn$lims, refn$lims, tolerance = tol,
             info = "bcanon: lims match baseline")
expect_equal(curn$stats, refn$stats, tolerance = tol,
             info = "bcanon: stats match baseline")
expect_equal(curn$equiv, refn$equiv, tolerance = tol,
             info = "bcanon: equiv match baseline")

## ---------------------------------------------------------------
## bcapar
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
curp <- bcapar(t0 = t0, tt = tt, bb = b_star)

expect_equal(curp$lims, refp$lims, tolerance = tol,
             info = "bcapar: lims match baseline")
expect_equal(curp$stats, refp$stats, tolerance = tol,
             info = "bcapar: stats match baseline")
expect_equal(curp$ustats, refp$ustats, tolerance = tol,
             info = "bcapar: ustats match baseline")
