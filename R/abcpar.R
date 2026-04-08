## Approximate Bootstrap Confidence (ABC) method for exponential families.
## Computes analytical BCa-like limits using numerical derivatives of func
## along eigenvectors of the sufficient statistic covariance.
## DiCiccio & Efron (1992) "More accurate confidence intervals in
## exponential families." Biometrika.
abcpar <- function (func,bb,lambda = 0.001,alpha = c(0.025, 0.05, 0.1, 0.16), ...)  {

    call <- match.call()

    b0 <- colMeans(bb)
    BB <- scale(bb, center = TRUE, scale = FALSE)
    y <- b0
    ## etahat: MLE of natural parameter; zero because bb is centered.
    ## Used as base point for exponential tilting in Mu().
    etahat <- numeric(length(b0))
    ## Exponential family mean function: tilts empirical distribution by exp(a'x).
    ## Maps natural parameter perturbation to mean parameter space. DE92 Sec 2.
    Mu <- function(a) {
        w <- as.vector(BB %*% a)
        W <- exp(w)
        W <- W / sum(W)
        as.vector(y + t(W %*% BB))}

    p <- length(y)
    S <- stats::var(bb)
    S <- as.matrix(S)
    thetahat <- func(y, ...)
    ## Eigen decomposition of sufficient statistic covariance.
    ## Derivatives are computed along eigenvectors (principal directions).
    eig <- eigen(S)
    if (any(eig$values <= 0))
        cli::cli_abort("Covariance matrix of {.arg bb} is singular or near-singular; ABC method requires full rank.")
    evals <- (eig$values)^0.5
    evecs <- as.matrix(eig$vectors)
    ## Numerical first (b.) and second (b..) derivatives of func along each
    ## eigendirection, scaled by sqrt(eigenvalue). DE92 Sec 3.
    b.. <- b. <- numeric(p)
    for (j in seq_len(p)) {
        b1 <- func(y + lambda * evals[j] * evecs[, j], ...)
        b2 <- func(y - lambda * evals[j] * evecs[, j], ...)
        b.[j] <- (b1 - b2) / (2 * lambda)
        b..[j] <- (b1 - 2 * thetahat + b2) / lambda^2
    }
    ## bhat = half sum of second derivatives (bias term)
    bhat <- sum(b..) / 2
    ## ehat = gradient in natural parameter space (scaled by inverse sqrt eigenvalues)
    ehat <- as.vector(evecs %*% (b./evals))
    ## dhat = gradient in mean parameter space (scaled by sqrt eigenvalues)
    dhat <- as.vector(evecs %*% (b. * evals))
    ## sighat = delta-method SE of theta-hat along focal direction
    sighat <- (sqrt(ehat %*% S %*% ehat))[1, 1]
    ## lam = step size normalized by SE for stable numerical differentiation
    lam <- lambda / sighat
    ## Acceleration via numerical 2nd derivative of Mu along ehat direction.
    ## The 6*sighat^3 denominator corresponds to (1/6) in the standard
    ## skewness-based acceleration formula. DE92 Sec 3.
    a0 <- sum(ehat * Mu(etahat))
    a1 <- sum(ehat * Mu(etahat + lam * ehat))
    a2 <- sum(ehat * Mu(etahat - lam * ehat))
    a <- (a1 - 2 * a0 + a2)/(lam^2 * 6 * sighat^3)
    ## delta = unit vector in focal direction (dhat normalized by SE)
    delta <- dhat / sighat
    ## cq = quadratic curvature of func along focal direction (finite difference).
    ## DE92 Sec 3.
    cq <- (func(y + lambda * delta, ...) + func(y - lambda * delta, ...) - 2 * thetahat) /
        (2 * sighat * lambda^2)
    ## curv = relative curvature: bias/SE minus quadratic curvature
    curv <- bhat / sighat - cq
    ## ABC bias-correction: acceleration minus curvature.
    ## Distinct from bootstrap-based z0 (which is empirical). DE92 Sec 3.
    z0 <- a - curv
    al <- c(alpha, .5, rev(1 - alpha))
    za <- stats::qnorm(al)
    ## BCa-like transformation with acceleration and bias correction
    z0a <- (z0 + za) / (1 - a * (z0 + za))
    ## Second-order correction including curvature
    z1a <- z0a + a * (z0a + z0)^2
    ## Standard (normal-theory) limits for comparison
    standard <- thetahat + sighat * za
    ## ABC limits: evaluate func at corrected quantile points
    ABC <- numeric(length(za))
    for (j in seq_along(za)) ABC[j] <- func(y + delta * z1a[j], ...)
    ## Quadratic approximation to ABC limits
    ABCquad <- thetahat + sighat * (z1a + cq * z1a^2)
    lims <- cbind( ABC, ABCquad, standard)
    dimnames(lims) <- list(al, c( "ABC", "ABCquad", "Stand"))
    list(call = call, lims = lims,
         thet.sig.bias = c(thetahat, sighat, bhat),
         a.z0.cq = c(a, z0, cq))
}
