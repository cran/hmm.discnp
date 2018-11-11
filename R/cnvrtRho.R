cnvrtRho <- function(Rho) {
    if(is.data.frame(Rho)) {
        ok <- ncol(Rho)==3 & identical(names(Rho)[1:2],c("y","state"))
        if(!ok)
            stop("Argument \"Rho\" is not of the correct form.\n")
        ly  <- levels(Rho$y)
        m   <- length(ly)
        n   <- length(levels(Rho$state))
        Rho <- matrix(exp(Rho$Intercept),nrow=m,ncol=n)
        Rho <- t(t(Rho)/apply(Rho,2,sum))
        rownames(Rho) <- ly
        return(Rho)
    }
    if(is.matrix(Rho)) {
        xxx <- vector("list",ncol(Rho))
        for(k in 1:ncol(Rho)) {
            xxx[[k]] <- p2expForm(Rho[,k])
        }
        rnms <- rownames(Rho)
        if(is.null(rnms)) rnms <- 1:nrow(Rho)
        y <- factor(rep(rnms,ncol(Rho)),levels=rnms)
        state <- factor(rep(1:ncol(Rho),each=nrow(Rho)))
        Intercept = do.call(c,xxx)
        return(data.frame(y=y,state=state,Intercept=Intercept))
    }
    stop("Argument \"Rho\" is not of an appropriate class.\n")
}    
