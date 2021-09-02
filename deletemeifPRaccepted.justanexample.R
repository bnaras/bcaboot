data(diabetes, package = "bcaboot")
Xy <- cbind(diabetes$x, diabetes$y)
rfun <- function(Xy) {
  y <- Xy[, 11]
  X <- Xy[, 1:10]
  summary(lm(y~X) )$adj.r.squared
}
set.seed(1234)
inds = matrix(NA, nrow = 1000, ncol = nrow(Xy))
inds = t(apply(inds, MARGIN = 1, function(inds) sample(1:nrow(Xy), replace = TRUE))) # get bootstrap indices to feed into estimator rfun()
tt = vector(mode = "numeric", length = 1000)
tt = apply(inds, MARGIN = 1, function(inds) rfun(Xy[inds,])) # get \hat{theta}*
Y = t(apply(inds, MARGIN = 1, tabulate, nbins = nrow(Xy)))
t0 = rfun(Xy)
bcaboot::bcajack2(B = list(Y=Y, tt=tt, t0=t0), alpha = c(0.025,0.975),
                  pct = 0.33)$lims[,1] # this produces NA BCa CIs, because a is NA, because B*pct=1000*0.33=330<442=n.
