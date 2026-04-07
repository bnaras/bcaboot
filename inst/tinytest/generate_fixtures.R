## Script to generate baseline regression fixtures.
## Run this ONCE on the current (pre-cleanup) code from the package root.
## The saved .rds files serve as the regression baseline.
##
## Usage: Rscript inst/tinytest/generate_fixtures.R

library(bcaboot)

fixture_dir <- "inst/tinytest/fixtures"
if (!dir.exists(fixture_dir)) dir.create(fixture_dir, recursive = TRUE)

## ---------------------------------------------------------------
## Common setup: diabetes data and rfun (used by bcajack, bcajack2, bcapar)
## ---------------------------------------------------------------
data(diabetes, package = "bcaboot")
Xy <- cbind(diabetes$x, diabetes$y)
rfun <- function(Xy) {
    y <- Xy[, 11]
    X <- Xy[, 1:10]
    summary(lm(y ~ X))$adj.r.squared
}

## --- bcajack ---
set.seed(1234)
ref_bcajack <- bcajack(x = Xy, B = 1000, func = rfun, m = 34, verbose = FALSE)
saveRDS(ref_bcajack, file.path(fixture_dir, "ref_bcajack.rds"))
cat("Saved ref_bcajack\n")

## --- bcajack2 ---
set.seed(1234)
ref_bcajack2 <- bcajack2(x = Xy, B = 1000, func = rfun, m = 40, verbose = FALSE)
saveRDS(ref_bcajack2, file.path(fixture_dir, "ref_bcajack2.rds"))
cat("Saved ref_bcajack2\n")

## --- bcanon ---
## Uses the logistic regression example (bcanon's smooth.spline diagnostic
## fails on the diabetes/rfun data — a pre-existing issue).
d <- readRDS("inst/examples/bcanon.RDS")
func5 <- function(Xy) {
    y <- Xy[, 7]; X <- Xy[, 1:6]
    re <- glm(y ~ X, binomial)
    re$coef[4]
}
set.seed(1234)
ref_bcanon <- bcanon(B = 1000, x = d$Xy, func = func5, verbose = FALSE)
saveRDS(ref_bcanon, file.path(fixture_dir, "ref_bcanon.rds"))
cat("Saved ref_bcanon\n")

## --- bcapar ---
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
set.seed(5678)
ref_bcapar <- bcapar(t0 = t0, tt = tt, bb = b_star)
saveRDS(ref_bcapar, file.path(fixture_dir, "ref_bcapar.rds"))
cat("Saved ref_bcapar\n")

cat("All fixtures saved to:", fixture_dir, "\n")
