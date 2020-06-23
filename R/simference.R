simference <- function(object,expForm=TRUE,seed=NULL,
                       nsim=100,verbose=TRUE) {
#
# Simulate multiple parameter estimates to use for Monte
# Carlo inference in respect of the parameter values.
#
    stopifnot(inherits(object,"hmm.discnp"))
    if(is.null(seed)) seed <- sample(1:1e5,1)
    set.seed(seed)
    xxx <- rhmm(object,nsim=nsim,verbose=verbose)
    if(verbose) {
        rslt <- vector("list",nsim)
        for(i in 1:nsim) {
            fit <- update(object,data=xxx[[i]])
            rslt[[i]] <- reparam(fit,expForm=expForm)
            cat(i,"")
            if(i%%10 == 0) cat("\n")
        }
        if(i%%10 != 0) cat("\n")
    } else {
        rslt <- lapply(xxx,function(x,obj) {
                               fit <- update(obj,data=x)
                               reparam(fit,expForm=expForm)
                           },obj=object)
    }
    M <- matrix(unlist(rslt),byrow=TRUE,nrow=nsim)
    colnames(M) = names(rslt[[1]])
    M
}
