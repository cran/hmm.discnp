cnvrtRho <- function(Rho) {
    if(inherits(Rho,"data.frame")) {
        ok <- ncol(Rho)==3 & identical(names(Rho),c("y","state","Intercept"))
        if(!ok)
            stop("Argument \"Rho\" is not of the correct form.\n")
        ly  <- levels(Rho$y)
        ls  <- levels(Rho$state)
        m   <- length(ly)
        n   <- length(ls)
        Rho <- matrix(Rho$Intercept,nrow=m,ncol=n)
        Rho <- apply(Rho,2,expForm2p)
        rownames(Rho) <- ly
        colnames(Rho) <- ls
        return(Rho)
    }
    if(inherits(Rho,"matrix")) {
        xxx <- vector("list",ncol(Rho))
        for(k in 1:ncol(Rho)) {
            xxx[[k]] <- p2expForm(Rho[,k])
        }
        rnms <- rownames(Rho)
        if(is.null(rnms)) rnms <- 1:nrow(Rho)
        y <- factor(rep(rnms,ncol(Rho)),levels=rnms)
        cnms <- colnames(Rho)
        if(is.null(cnms)) cnms <- 1:ncol(Rho)
        state <- factor(rep(cnms,each=nrow(Rho)))
        Intercept = do.call(c,xxx)
        return(data.frame(y=y,state=state,Intercept=Intercept))
    }
    stop("Argument \"Rho\" is not of an appropriate class.\n")
}    
