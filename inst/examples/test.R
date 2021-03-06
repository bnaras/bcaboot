dir <- "../../R/"
source(paste0(dir, "bcajack.R"))
source(paste0(dir, "bcajack2.R"))
source(paste0(dir, "bcapar.R"))
source(paste0(dir, "bcaplot.R"))

library(bcaboot)

library(lars)
data(diabetes)
Xy <- cbind(diabetes$x, diabetes$y)
rfun <- function(Xy) {
    y <- Xy[, 11]
    X <- Xy[, 1:10]
    summary(lm(y~X) )$adj.r.squared
}
m <- bcajack(x = Xy, B = 2000, func = rfun)

bcajack(x = Xy, B = 2000, func = rfun, m = 40)
###
###
library(bcaboot)
bcajack(x = Xy, B = 2000, func = rfun, m = 40, catj = 0)
library(lars)
data(diabetes)
Xy <- cbind(diabetes$x, diabetes$y)
rfun <- function(Xy) {
    y <- Xy[, 11]
    X <- Xy[, 1:10]
    summary(lm(y~X) )$adj.r.squared
}
bcajack(x = Xy, B = 2000, func = rfun)

bcajack(x = Xy, B = 2000, func = rfun, m = 40)

bcajack2(x = Xy, B = 2000, func = rfun)

r1 <- bcajack(x = Xy, B = 2000, func = rfun)
