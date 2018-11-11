squantCI <- function(object,seed=NULL,alpha=0.05,nsim=100,verbose=TRUE) {
#
# Estimated confidence intervals for model parmeters, based on quantiles
# of estimates obtained via simulation.
#
    stopifnot(inherits(object,"hmm.discnp"))
    if(alpha <= 0 | alpha >= 1)
        stop("Argument \"alpha\" must be strictly between 0 and 1.\n")
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
    CIs <- apply(M,2,function(x,a){quantile(x,c(a/2,1-a/2))},a=alpha)
    attr(CIs,"seed") <- seed
    return(CIs)
}
