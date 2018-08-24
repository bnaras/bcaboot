library(bcaboot)

set.seed(3891)
B <- 16000

sigma1 <- sqrt(2)
sigma2 <- sqrt(2)
n1 <- 10
n2 <- 42

sample1 <- sigma1^2 * rchisq(n = 100, df = n1) / n1
sample2 <- sigma2^2 * rchisq(n = 100, df = n2) / n2

v1 <- var(sample1)
v2 <- var(sample2)

do_ratio_boot <- function(B, v1, v2) {
    s1 <- sqrt(v1) * rchisq(n = B, df = n1)  / n1
    s2 <- sqrt(v2) * rchisq(n = B, df = n2)  / n2
    theta_star <- s1 / s2
    beta_star <- cbind(s1, s2)
    list(theta = v1 / v2,
         theta_star = theta_star,
         bb = beta_star)
}

ratio_stuff <- do_ratio_boot(B, 1, 1)
funcF <- function(beta) {
    beta[1] / beta[2]
}
result <- bcapar(t0 = ratio_stuff$theta,
                 tt = ratio_stuff$theta_star,
                 bb = ratio_stuff$bb, f = funcF)
