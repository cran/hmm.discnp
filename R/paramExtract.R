paramExtract <- function(Rho) {
    lvls <- levels(Rho$y)
    m    <- length(lvls)
    lm   <- lvls[m]
    ok   <- Rho$y!=lm
    as.vector(as.matrix(Rho[ok,-(1:2)]))
}
