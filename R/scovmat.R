scovmat <- function(object,expForm=TRUE,seed=NULL,
                    nsim=100,verbose=TRUE) {
#
# Simulated covariance matrix; i.e. an estimate of the covariance
# matrix of the parmeter estimates obtained via simulation.
#
    M   <- simference(object=object,expForm=expForm,seed=seed,
                      nsim=nsim,verbose=verbose)
    ccc <- var(M)
    attr(ccc,"seed") <- attr(M,"seed")
    ccc
}
