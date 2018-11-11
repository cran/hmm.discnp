scovmat <- function(object,seed=NULL,nsim=100,verbose=TRUE) {
#
# Simulated covariance matrix; i.e. an estimate of the covariance
# matrix of the parmeter estimates obtained via simulation.
#
    stopifnot(inherits(object,"hmm.discnp"))
    if(is.null(seed)) seed <- sample(1:1e5,1)
    set.seed(seed)
    newstyle <- object$newstyle
    ylengths <- object$ylengths
    xxx <- rhmm(object,ylengths=ylengths,nsim=nsim,verbose=verbose)
    if(verbose) {
        rslt <- vector("list",nsim)
        for(i in 1:nsim) {
            fit <- update(object,data=xxx[[i]])
            if(newstyle) fit$Rho <- cnvrtRho(fit$Rho)
            rslt[[i]] <- reparam(fit)
            cat(i,"")
            if(i%%10 == 0) cat("\n")
        }
        if(i%%10 != 0) cat("\n")
    } else {
        rslt <- lapply(xxx,function(x,obj,newstyle) {
                               fit <- update(obj,data=x)
                               if(newstyle) fit$Rho <- cnvrtRho(fit$Rho)
                               reparam(fit)
                           },obj=object,newstyle=newstyle)
    }
    M <- matrix(unlist(rslt),byrow=TRUE,nrow=nsim)
    colnames(M) = names(rslt[[1]])
    ccc <- var(M)
    attr(ccc,"seed") <- seed
    return(ccc)
}
