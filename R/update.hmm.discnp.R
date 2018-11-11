update.hmm.discnp <- function(object,...,data,Kplus1=FALSE,tpm2=NULL,
                              verbose=FALSE,method=NULL,optimiser=NULL,
                              stationary=NULL,mixture=NULL,cis=NULL,
                              tolerance=NULL,itmax=NULL,crit=NULL,X=NULL,
                              addIntercept=NULL) {

# Check that the "data" argument is compatible with the "object"
# argument.
checkDat  <- tidyList(data)
if(attr(checkDat,"parity") != object$parity)
    stop("The \"object\" and \"data\" arguments are not compatible.\n")

if(is.null(method))       method       <- object$args$method
if(is.null(optimiser))    optimiser    <- object$args$optimiser
if(is.null(stationary))   stationary   <- object$args$stationary
if(is.null(mixture))      mixture      <- object$args$mixture
if(is.null(cis))          cis          <- object$args$cis
if(is.null(tolerance))    tolerance    <- object$args$tolerance
if(is.null(itmax))        itmax        <- object$args$itmax
if(is.null(crit))         crit         <- object$args$crit
if(is.null(addIntercept)) addIntercept <- object$args$addIntercept

# Increment the number of states by one if required.
if(Kplus1) {
    tpm <- object$tpm
    if(identical(tpm,NA)) {
        K   <- 1
        tpm <- if(is.null(tpm2)) matrix(0.5,2,2) else tpm2
        if(!identical(dim(tpm),c(2L,2L))) {
            stop("Argument \"tpm2\" must be a 2 x 2 matrix.\n")
        }
    } else {
        K   <- nrow(tpm)
        tpm[,K] <- tpm[,K]/2
        tpm <- cbind(tpm,tpm[,K])
        tpm <- rbind(tpm,tpm[K,])
    }
    if(cis) {
        ispd <- object$ispd
        if(identical(ispd,NA)) {
            ispd <- c(0.5,0.5)
        } else if(stationary) {
            ispd <- revise.ispd(tpm)
        } else {
            K <- length(ispd)
            ispd[K] <- ispd[K]/2
            ispd <- c(ispd,ispd[K])
        }
    } else ispd <- NULL
    Rho  <- object[["Rho"]]
    nval <- length(levels(Rho[["y"]]))
    xxx  <- Rho[(1+(K-1)*nval):(K*nval),]
    xxx$state <- factor(K+1)
    Rho  <- rbind(Rho,xxx)
    par0 <- list(ispd=ispd,tpm=tpm,Rho=Rho)
} else {
    par0 <- object[c("ispd","tpm","Rho")]
    Rho  <- object[["Rho"]]
}

# Deal with the predictor variables if these are supplied.
if(object$args$newstyle && ncol(Rho) > 3) {
# X used previously:
    if(is.null(X)) stop("Argument \"X\" must be specified.\n")
    X <- tidyList(X,rp="predictor",addIntercept=addIntercept)
    prednames <- attr(X,"prednames")
    OK <- identical(prednames,colnames(Rho)[-(1:2)])
    if(!OK) {
        nnn <- paste(prednames,collapse=", ")
        ttt <- paste(colnames(Rho)[-(1:2)],collapse=", ")
        whinge <- paste0("The \"X\" that was supplied yields predictor",
                         " names:\n",nnn,"\n",
                         "The previous predictor names were:\n",
                         ttt,"\n",
                         "What we've got here is a failure to be",
                         " compatible.\n")
        stop(whinge)
    }
} else {
# X being added to the mix:
    if(!is.null(X)) {
        X <- tidyList(X,rp="predictor",addIntercept=addIntercept)
        prednames <- attr(X,"prednames")
        np        <- length(prednames)
        newRho    <- cbind(Rho,matrix(0,nrow=nrow(Rho),ncol=np-1))
        names(newRho)[3:(np+2)] <- prednames
        par0[["Rho"]] <- newRho
    }
}

ill <- logLikHmm(y=data,X=X,tpm=par0$tpm,Rho=par0$Rho)

if(identical(par0$tpm,NA)) {
    K <- 1
    par0 <- NULL
} else K <- NULL

fit <- hmm(y=data,K=K,par0=par0,method=method,optimiser=optimiser,
           stationary=stationary,mixture=mixture,cis=cis,
           tolerance=tolerance,itmax=itmax,crit=crit,X=X,
           addIntercept=addIntercept,verbose=verbose,...)
nafter <- which(names(fit)=="log.like") - 1
fit <- append(fit,list(init.log.like=ill),after=nafter)
class(fit) <- "hmm.discnp"
return(fit)
}
