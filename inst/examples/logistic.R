library(bcaboot)
library(glmnet)
load("neonates.rda")

X <- as.matrix(neonates[, seq_len(11)])
Y <- neonates$y
weights <- ifelse(Y == 0, 1, 4)
neo_glm <- glm(formula = y ~ ., family = "binomial", weights = weights, data = neonates)

set.seed(3891)
B <- 2000

fit_model <- function(data, var = "resp") {
    ## response should be named y
    coef(glm(formula = y ~ ., data = data, family = "binomial"))[var]
}

do_boot <- function(B, glm_model, var = "resp") {
    pi_hat <- glm_model$fitted.values
    n <- length(pi_hat)
    y_star <- sapply(seq_len(B), function(i) ifelse(runif(n) <= pi_hat, 1, 0))
    beta_star <- apply(y_star, 2,
                       function(y) {
                           boot_data <- glm_model$data
                           boot_data$y <- y
                           coef(glm(formula = y ~ ., data = boot_data, weights = weights, family = "binomial"))
                       })
    list(theta = coef(glm_model)[var],
         theta_star = beta_star[var, ],
         bb = t(y_star) %*% model.matrix(glm_model),
         y_star = y_star)
}

glm_stuff <- do_boot(B, neo_glm)
glm_bca <- bcapar(t0 = glm_stuff$theta,
                  tt = glm_stuff$theta_star,
                  bb = glm_stuff$bb)


set.seed(3891)
B <- 2000

neo_glmnet <- cv.glmnet(x = X, y = Y, family = "binomial", weights = weights)
coef <- as.matrix(coef(neo_glmnet, s = neo_glmnet$lambda.min))

do_glmnet_boot <- function(B, X, y, glmnet_model, var = "resp") {
    lambda <- glmnet_model$lambda.min
    ## preds <- predict(glmnet_model, newx = X, s = "lambda.min")
    ## pi_hat <- exp(preds) / (1 + exp(preds))
    pi_hat <- predict(glmnet_model, newx = X, s = "lambda.min", type = "response")
    n <- length(pi_hat)
    y_star <- sapply(seq_len(B), function(i) ifelse(runif(n) <= pi_hat, 1, 0))
    beta_star <- apply(y_star, 2,
                       function(y) {
                           as.matrix(coef(glmnet(x = X, y = y, lambda = lambda, weights = weights, family = "binomial")))
                       })
    coefs <- as.matrix(coef(glmnet_model, s = glmnet_model$lambda.min))
    rownames(beta_star) <- rownames(coefs)
    list(theta = coefs[var, ],
         theta_star = beta_star[var, ],
         bb = t(y_star) %*% X,
         y_star = y_star)
}

glmnet_stuff <- do_glmnet_boot(B, X, y, neo_glmnet)

glmnet_bca <- bcapar(t0 = glmnet_stuff$theta,
                     tt = glmnet_stuff$theta_star,
                     bb = glmnet_stuff$bb)
bcaplot(glmnet_bca)

