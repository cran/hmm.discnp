squantCI <- function(object,expForm=TRUE,seed=NULL,alpha=0.05,
                     nsim=100,verbose=TRUE) {
#
# Estimated confidence intervals for model parmeters, based on quantiles
# of estimates obtained via simulation.
#
    M   <- simference(object=object,expForm=expForm,seed=seed,
                      nsim=nsim,verbose=verbose)
    CIs <- t(apply(M,2,function(x,a){quantile(x,c(a/2,1-a/2))},a=alpha))
    attr(CIs,"seed") <- seed
    return(CIs)
}
