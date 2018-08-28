abcpar <- function (func,bb,lambda = 0.001,alpha = c(0.025, 0.05, 0.1, 0.16), ...)  {

    call <- match.call()

    b0 <- colMeans(bb)
    BB <- scale(bb, center = TRUE, scale = FALSE)
    y <- b0
    etahat <- numeric(length(b0))
    Mu <- function(a) {
        w <- as.vector(BB %*% a)
        W <- exp(w)
        W <- W / sum(W)
        as.vector(y + t(W %*% BB))}

    p <- length(y)
    S <- stats::var(bb)
    S <- as.matrix(S)
    I <- diag(p)
    thetahat <- func(y, ...)
    eig <- eigen(S)
    evals <- (eig$values)^0.5
    evecs <- as.matrix(eig$vectors)
    b.. <- b. <- numeric(p)
    for (j in seq_len(p)) {
        b1 <- func(y + lambda * evals[j] * evecs[, j], ...)
        b2 <- func(y - lambda * evals[j] * evecs[, j], ...)
        b.[j] <- (b1 - b2) / (2 * lambda)
        b..[j] <- (b1 - 2 * thetahat + b2) / lambda^2
    }
    bhat <- sum(b..) / 2
    ehat <- as.vector(evecs %*% (b./evals))
    dhat <- as.vector(evecs %*% (b. * evals))
    ## sighat <- sqrt(ehat %*% S %*% ehat)
    sighat <- (sqrt(ehat %*% S %*% ehat))[1, 1]
    lam <- lambda / sighat
    a0 <- sum(ehat * Mu(etahat))
    a1 <- sum(ehat * Mu(etahat + lam * ehat))
    a2 <- sum(ehat * Mu(etahat - lam * ehat))
    a <- (a1 - 2 * a0 + a2)/(lam^2 * 6 * sighat^3)
    delta <- dhat / sighat
    cq <- (func(y + lambda * delta, ...) + func(y - lambda * delta, ...) - 2 * thetahat) /
        (2 * sighat * lambda^2)
    curv <- bhat / sighat - cq
    z0 <- a - curv
    al <- c(alpha, .5, rev(1 - alpha))
    za <- stats::qnorm(al)
    z0a <- (z0 + za) / (1 - a * (z0 + za))
    z1a <- z0a + a * (z0a + z0)^2
    standard <- thetahat + sighat * za
    ABC <- numeric(length(za))
    for (j in seq_along(za)) ABC[j] <- func(y + delta * z1a[j], ...)
    ABCquad <- thetahat + sighat * (z1a + cq * z1a^2)
    lims <- cbind( ABC, ABCquad, standard)
    dimnames(lims) <- list(al, c( "ABC", "ABCquad", "Stand"))
    list(call = call, lims = lims,
         thet.sig.bias = c(thetahat, sighat, bhat),
         a.z0.cq = c(a, z0, cq))
}
