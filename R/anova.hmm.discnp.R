anova.hmm.discnp <- function(object,...){
    obs  <- c(list(object),list(...))
    if(length(obs) != 2)
        stop("This anova method handles only comparisons between two models.\n")
    ok <- inherits(obs[[2]],"hmm.discnp")
    if(!all(ok)) stop("Second argument is not of class \"hmm.discnp\".\n")
    np <- sapply(obs,function(x){x$npar})
    o  <- order(np)
    obs <- obs[o]
    np <- np[o]
    nu   <- np[2]-np[1]
    ll <- sapply(obs,function(x){x$log.like})
    stat <- 2*(ll[2]-ll[1])
    pv   <- pchisq(stat,nu,lower.tail=FALSE)
    rslt <- list(stat=stat,df=nu,pvalue=pv)
    attr(rslt,"details") <- c(ll1=ll[1],ll2=ll[2],np1=np[1],np2=np[2])
    rslt
}
