hmm <- function(y,yval=NULL,par0=NULL,K=NULL,rand.start=NULL,
                method=c("EM","bf","LM","SD"),hglmethod=c("fortran","oraw","raw"),
                optimiser=c("nlm","optim"),optimMethod=NULL,stationary=cis,
                mixture=FALSE,cis=TRUE,indep=NULL,tolerance=1e-4,digits=NULL,
                verbose=FALSE,itmax=200,crit=c("PCLL","L2","Linf","ABSGRD"),
                keep.y=FALSE,keep.X=keep.y,newstyle=TRUE,X=NULL,
                addIntercept=TRUE,lmc=10,hessian=FALSE,...) {

#
# Function hmm.  To fit a Hidden Markov model to a data set where the
# observations come from one of a number of finite discrete
# distributions, depending on the (hidden) state of the Markov chain.
# These distributions are specified by a matrix Rho = [rho_ij] where
# rho_ij = P(Y = y_i | S = j), Y being the observable random variable
# and S being the hidden state.
# 
# Note that y used to be allowed to be a matrix, each column
# being interpreted as an independent replicate of the observation
# sequence.  Now y should be either a vector or a one or two column
# matrix or a *list* of vectors or a *list* of one or two column
# matrices.
#

method <- match.arg(method)
# EM <--> Expectation/maximisation algorithm.
# bf <--> brute force (maximise using either nlm() or optim()).
# LM <--> Levenberg-Marquardt algorithm.
# SD <--> steepest descent

# Check the "crit" argument.
crit <- match.arg(crit)
if(method=="EM" & crit=="ABSGRD")
    stop("Stopping criterion \"ABSGRD\" is not useable with the EM algorithm.\n")

if(!newstyle) {
    if(method %in% c("bf","LM","SD")) {
        whinge <- paste0("When the \"",method,"\" method is used ",
                         "\"newstyle\" must be TRUE.\n")
        stop(whinge)
    }
}

if(method %in% c("LM","SD")) {
    if(!is.null(X)){
        whinge <- paste0("Cannot currently use auxiliary predictors with the ",
                         method," method.\n")
        stop(whinge)
    }
    if(!stationary) {
        whinge <- paste0("Currently the ",method," method only handles",
                         " stationary models.\n")
        stop(whinge)
    }
}

# Check on consistency of ``mixture'' and ``stationary''.
# Also check on "mixture" and "method".
if(mixture) {
    if(!stationary)
	stop("Makes no sense for a mixture model to be non-stationary.\n")
    if(method=="bf")
        stop("Currently the \"bf\" method does not do mixtures.\n")
    if(method=="LM")
        stop("Currently the \"LM\" method does not do mixtures.\n")
}

# Check on consistancy of ``stationary'' and ``cis''.
if(stationary & !cis)
	stop(paste("If the model is stationary the initial state\n",
                   "probability distribution must be constant\n",
                   "across data sequences.\n"))

# Check that y is sensibly consistent; if it is a single
# vector or matrix, put it into a single-entry list.
# In some circumstances "y" may have already been processed
# into a "tidyList" object; the tidyList() function has now
# (27/18/2018) been adjusted so that if this is the case then "y"
# is returned unchanged.  This avoids (amongst other things?)  the
# problem that if tidyList() has been applied then the entries of
# "y" will have been coerced to character mode.  Consequently the
# "numeric" attribute of the result could be wrong if tidyList()
# were re-applied.
y <- tidyList(y)
bivar <- attr(y,"parity")=="bivar"
if(bivar) {
    if(method %in% c("bf","LM")) {
        whinge <- paste0("Currently the \"",method,"\" method can only ",
                         "be used for univariate data.\n")
        stop(whinge)
    }
}

# Check that at least one of par0 and K is specified and, if
# both are specified, there is no inconsistency.
if(is.null(K)) {
    if(is.null(par0))
	stop('One of par0 and K must be specified.')
    K <- nrow(par0$tpm)
    if(K != ncol(par0$tpm))
            stop("The specified tpm is not square.\n")
} else if(K==1 & !is.null(par0)) {
    warning("When K equals 1, par0 is ignored.\n")
} else if(!is.null(par0) && K != nrow(par0$tpm)) {
    stop("The values of \"K\" and \"par0\" are inconsistent.\n")
}

# Check on (conditional) independence requirement.
if(bivar) {
    if(is.null(par0$Rho)) {
        if(is.null(indep)) {
            stop("Neither \"indep\" nor \"par0$Rho\" have been supplied.\n")
        }
    } else if(is.list(par0$Rho)) {
        if(is.null(indep)) {
            indep <- TRUE
        } else if(!indep) {
             stop(paste("The value of \"indep\" and that of \"par0$Rho\"\n",
                        "are inconsistent.\n",sep=""))
        }
    } else {
        if(is.null(indep)) {
            indep <- FALSE
        } else if(indep) {
             stop(paste("The value of \"indep\" and that of \"par0$Rho\"\n",
                        "are inconsistent.\n",sep=""))
        }
    }
}

# Set up the multiplier for BIC.
nobs <- sum(!is.na(unlist(y)))
if(bivar) nobs <- nobs/2
bicm <- log(nobs)

# BI <--> "bivariate independent".
# BD <--> "bivariate dependent".
# UV <--> "univariate".
if(bivar) {
    if(indep) {
        rslt <- hmmBI(y,yval=yval,par0=par0,K=K,rand.start=rand.start,
                      stationary=stationary,mixture=mixture,cis=cis,
                      tolerance=tolerance,verbose=verbose,itmax=itmax,
                      crit=crit,bicm=bicm)
    } else {
        rslt <- hmmBD(y,yval=yval,par0=par0,K=K,rand.start=rand.start,
                      stationary=stationary,mixture=mixture,cis=cis,
                      tolerance=tolerance,verbose=verbose,itmax=itmax,
                      crit=crit,bicm=bicm)
    }
} else {
    optimiser <- match.arg(optimiser)
    rslt <- hmmUV(y,yval=yval,par0=par0,K=K,rand.start=rand.start,
                  method=method,hglmethod=hglmethod,optimiser=optimiser,
                  optimMethod=optimMethod,stationary=stationary,
                  mixture=mixture,cis=cis,tolerance=tolerance,
                  digits=digits,verbose=verbose,itmax=itmax,crit=crit,
                  bicm=bicm,newstyle=newstyle,X=X,
                  addIntercept=addIntercept,lmc=lmc,hessian=hessian,...)
}
ylnths <- sapply(y,nrow)
if("prior.emsteps" %in% names(rslt)) {
    naft <- which(names(rslt)=="prior.emsteps")
} else {
    naft <- which(names(rslt)=="nstep")
}
rslt   <- append(rslt,list(ylengths=ylnths),after=naft)
naft <- naft + 1
nafrac <- nafracCalc(y)
rslt   <- append(rslt,list(nafrac=nafrac),after=naft)
naft <- naft + 1
if(keep.y) {
    rslt   <- append(rslt,list(y=y),after=naft)
    naft <- naft + 1
}
keep.X <- keep.X & !bivar & newstyle & !is.null(X)
if(keep.X) {
    rslt   <- append(rslt,list(X=X),after=naft)
    naft <- naft + 1
}
rslt   <- append(rslt,list(parity=attr(y,"parity")),after=naft)
naft <- naft + 1
rslt   <- append(rslt,list(numeric=attr(y,"numeric")),after=naft)
args <- list(newstyle=newstyle,method=method,optimiser=optimiser,
             optimMethod=optimMethod,stationary=stationary,
             mixture=mixture,cis=cis,tolerance=tolerance,
             itmax=itmax,crit=crit,addIntercept=addIntercept)
rslt$args   <- args
class(rslt) <- "hmm.discnp"
rslt
}
