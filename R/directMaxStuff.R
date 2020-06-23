#
# Functions that are needed when method is direct maximisation
# (via nlm() or optim()).
#

p2expForm <- function (x) {
# Convert a vector of probabilities (summing to 1) to
# a vector of parameters for the logistic style
# parametrisation of these probabilities.  Notice that
# the last entry of this vector is constrained to be 0.
# There has to be a constraint of course, corresponding
# to the original constraint that the probabilities
# sum to 1.
    z  <- log(x)
    nok <- z==-Inf
    if(any(nok)) {
        btm <- min(z[!nok])
        z[nok] <- min(btm,-300)
    }
    z <- z - z[length(z)]
    return(z)
}

expForm2p <- function(x){
# Convert a vector of parameters for the logistic style
# parameterisation to a vector of probabilities summing
# to 1.
    m  <- max(x)
    xr <- exp(x-m)
    xr/sum(xr)
}

reparam <- function(object,expForm=TRUE,stationary=NULL) {
# Convert a list with entries "ispd", "tpm" and "Rho", and possibly
# "stationary" to a vector of parameters (with no redundancy and
# the same information content).  Note that "object" may be an
# object of class "hmm.discnp" as returned by the function of
# that name.

if(is.null(stationary)) stationary <- object$stationary
if(is.null(stationary))
    stop("No value of \"stationary\" provided.\n")

newpars <- vector("list",3)
# Do ispd if "stationary" is not TRUE.
if(stationary) {
   newpars[[1]] <- NULL
} else {
   ispd  <- object$ispd
   K     <- length(ispd)
   if(expForm) {
       xxx   <- p2expForm(ispd)[-K]
       names(xxx) <- paste0("omega.",1:K)[-K]
   } else {
       xxx <- ispd[-K]
       names(xxx) <- paste0("pi.",1:K)[-K]
   }
   newpars[[1]] <- xxx
}

# Do tpm:
    tpm <- object$tpm
    K   <- ncol(tpm)
if(expForm) {
    nms <- paste("zeta",row(tpm)[,-K],col(tpm)[,-K],sep=".")
    A   <- t(apply(tpm,1,p2expForm))[,-K]
} else {
    nms <- paste("p",row(tpm)[,-K],col(tpm)[,-K],sep=".")
    A   <- tpm[,-K]
}
# The matrix A, of "zeta"  values or of "p" values, is strung out,
# *column by column*.
    xxx <- as.vector(A)
    names(xxx) <- nms
    newpars[[2]] <- xxx

# Do Rho:
    Rho  <- object$Rho
    if(!inherits(Rho,"data.frame"))
        stop("Component \"Rho\" of \"object\" is not of the correct form.\n")
    lvls <- levels(Rho$y)
    m    <- length(lvls)
    sts  <- unique(Rho$state)

# Changed 17/03/2018 to treat the *last* row of coefficients,
# for all states (last phi value, phi_{m,j}, all j) as 0.
# Rather than the first row.
if(expForm) {
    Rho  <- Rho[Rho$y!=lvls[m],-(1:2),drop=FALSE]
    cns  <- colnames(Rho)
    if(length(cns)==1) cns <- "phi" # Not "Intercept"!
    eee  <- expand.grid(lvls[-m],sts,cns)
    nms  <- do.call(paste,c(eee[,c(3,1,2)],list(sep=".")))
    Rho  <- as.matrix(Rho)
} else {
    if(ncol(Rho) > 3)
        stop("When there are predictors, only expForm=TRUE is allowed.\n")
    Rho <- cnvrtRho(Rho)
    Rho <- Rho[-m,]
    eee <- as.matrix(expand.grid(rownames(Rho)[-K],colnames(Rho)))
    nms <- apply(eee,1,function(x){paste0("rho.",paste(x,collapse="."))})
}

# The matrix of coefficient values is strung out *column by column*.
xxx  <- as.vector(Rho)
names(xxx) <- nms
newpars[[3]] <- xxx

# Return the result.
unlist(newpars)
}

getIspd <- function(pars,K) {
    return(expForm2p(c(pars[1:(K-1)],0)))
}

getTpm <- function(pars,K,stationary) {
# Get the transition probability matrix from the
# vector of (non-redundant) vector of parameters
# of the model.
    tpmind <- 1:(K*(K-1)) + if(stationary) 0 else K-1
    zeta <- cbind(matrix(pars[tpmind],nrow=K),0)
    tpm  <- t(apply(zeta,1,expForm2p))
    return(tpm)
}

getRho <- function(pars,K,rhovals,stationary,prednames) {
    lout   <- (K-1)*(if(stationary) K else K+1)
    outind <- 1:lout
    m      <- length(rhovals)
    rose   <- matrix(pars[-outind],nrow=K*(m-1))
    if(ncol(rose) != length(prednames)) {
        stop("Wrong length for \"pars\".\n")
    }
    y      <- factor(rep(rhovals,K),levels=rhovals)
    M      <- matrix(0,nrow=K*m,ncol=ncol(rose))
    M[y!=rhovals[m],] <- rose
    Rho  <- data.frame(y=y,state=factor(rep(1:K,each=length(rhovals))),M)
    colnames(Rho)[3:ncol(Rho)] <- prednames
    return(Rho)
}
