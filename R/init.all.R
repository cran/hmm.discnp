init.all <- function(nval,K,rand.start,mixture,indep,newstyle,yval,
                     prednames) {
#
# Function init.all to create (rather arbitrary) initial values for
# tpm and Rho.
#

# Do tpm --- the easy bit!!!
if(rand.start$tpm)
	tpm <- if(mixture) matrix(runif(K),K,K,byrow=TRUE)
			else matrix(runif(K*K),K,K)
else tpm <- matrix(1/K,K,K) + 1.5*diag(K)
tpm <- tpm/apply(tpm,1,sum)

# Now do Rho --- the hard bit.
if(length(nval)==1) {
    hasInt <- prednames[1] == "Intercept"
    ncX <- length(prednames)
    if(hasInt) ncX <- ncX - 1
    nyv <- length(yval)
    y   <- factor(rep(yval,K),levels=yval)
    ss  <- factor(rep(1:K,each=nyv),levels=1:K)
    if(ncX) {
        if(rand.start$Rho) {
            xxx <- vector("list",K)
            for(k in 1:K) {
                xxx[[k]] <- rbind(matrix(rnorm((nyv-1)*ncX),ncol=ncX),0)
            }
            cX <- do.call(rbind,xxx)
        } else cX <- matrix(0,nrow=nyv*K,ncol=ncX)
    } else cX <- NULL
    if(hasInt) {
        ccc <- vector("list",K)
        for(k in 1:K) {
            ccc[[k]] <- if(rand.start$Rho) c(rnorm(nyv-1),0) else
                            p2expForm(seq(from=k,by=K,length=nyv))
        }
        beta0 <- do.call(c,ccc)
    } else beta0 <- NULL
    cX           <- cbind(beta0,cX)
    colnames(cX) <- prednames
    Rho          <- data.frame(y=y,state=ss,cX)
    if(!newstyle) Rho <- cnvrtRho(Rho)
    return(list(tpm=tpm,Rho=Rho))
} else {
    if(indep) {
        Rho <- vector("list",2)
        for(i in 1:2) {
            if(rand.start$Rho) Rho[[i]] <- matrix(runif(K*nval[i]),K,nval[i])
            else Rho[[i]] <- matrix(1:(nval[i]*K),K,nval[i])
            Rho[[i]] <- t(Rho[[i]]/apply(Rho[[i]],1,sum))
        }
    } else {
        if(rand.start$Rho)
             Rho <- array(runif(K*prod(nval)),dim=c(nval[1],nval[2],K))
        else
             Rho <- array(1:(prod(nval)*K),dim=c(nval[1],nval[2],K))
        div <- apply(Rho,3,sum)
        Rho <- aperm(Rho,c(3,1,2))/div
        Rho <- aperm(Rho,c(2,3,1))
    }
}

list(tpm=tpm,Rho=Rho)
}
