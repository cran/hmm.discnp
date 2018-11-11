paramExtract <- function(Rho,newstyle) {
    if(newstyle) {
        lvls <- levels(Rho$y)
        m    <- length(lvls)
        lm   <- lvls[m]
        ok   <- Rho$y!=lm
        return(as.vector(as.matrix(Rho[ok,-(1:2)])))
    }
# Here Rho is (always) a matrix, so as.matrix() is not needed.
    as.vector(Rho[-nrow(Rho),])
}
