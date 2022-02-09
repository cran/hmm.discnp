hmm <- function(y,yval=NULL,par0=NULL,K=NULL,rand.start=NULL,
                method=c("EM","bf","LM","SD"),hglmethod=c("fortran","oraw","raw"),
                optimiser=c("nlm","optim"),optimMethod=NULL,stationary=cis,
                mixture=FALSE,cis=TRUE,indep=NULL,tolerance=1e-4,digits=NULL,
                verbose=FALSE,itmax=200,crit=c("PCLL","L2","Linf","ABSGRD"),
                X=NULL,keep.y=FALSE,keep.X=keep.y,
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

# If y is of class "multipleHmmDataSets" bail out.
if(inherits(y,"multipleHmmDataSets")) {
    whinge <- paste0("Argument \"y\" is a list of \"appropriate\" data sets.\n",
                     "  This function should be applied to the individual\n",
                     "  components, perhaps by means of lapply().\n")
    stop(whinge)
}

# EM <--> Expectation/maximisation algorithm.
# bf <--> brute force (maximise using either nlm() or optim()).
# LM <--> Levenberg-Marquardt algorithm.
# SD <--> steepest descent

# Check the "crit" argument.
crit <- match.arg(crit)
if(method=="EM" & crit=="ABSGRD")
    stop("Stopping criterion \"ABSGRD\" is not useable with the EM algorithm.\n")

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

# Check on consistency of "mixture" with "stationary" and "method".
if(mixture) {
    if(!stationary)
	stop("Makes no sense for a mixture model to be non-stationary.\n")
    if(method=="bf")
        stop("Currently the \"bf\" method does not do mixtures.\n")
    if(method=="LM")
        stop("Currently the \"LM\" method does not do mixtures.\n")
}

# Check on consistancy of "stationary" and "cis".
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
y     <- tidyList(y,yval=yval)
lvls  <- attr(y,"lvls")

# Check on method and predictors w.r.t. 
bivar <- attr(y,"parity")=="bivar"
if(bivar) {
    if(method %in% c("bf","LM")) {
        whinge <- paste0("Currently the \"",method,"\" method can only ",
                         "be used for univariate data.\n")
        stop(whinge)
    }
} else {
    if(!is.null(X) & method %in% c("EM","bf")) {
        X <- tidyList(X,rp="predictor",addIntercept=addIntercept)
        checkyXoK(y,X)
    }
    if(!is.null(X)) {
        prednames <- attr(X,"prednames")
    } else {
        prednames <- "Intercept"
    }
}

# Organise the starting values.
starts <- checkStartVal(par0,K,indep,lvls,rand.start,
                      mixture,prednames)
par0   <- starts$par0
K      <- starts$K
indep  <- starts$indep
stnms  <- starts$stnms

# Set up the multiplier for BIC.
nobs <- sum(!is.na(unlist(y)))
if(bivar) nobs <- nobs/2
bicm <- log(nobs)

# BI <--> "bivariate independent".
# BD <--> "bivariate dependent".
# UV <--> "univariate".
if(bivar) {
    if(indep) {
        rslt <- hmmBI(y,par0=par0,K=K,
                      stationary=stationary,mixture=mixture,cis=cis,
                      tolerance=tolerance,digits=digits,verbose=verbose,itmax=itmax,
                      crit=crit,bicm=bicm)
    } else {
        rslt <- hmmBD(y,par0=par0,K=K,
                      stationary=stationary,mixture=mixture,cis=cis,
                      tolerance=tolerance,digits=digits,verbose=verbose,itmax=itmax,
                      crit=crit,bicm=bicm)
    }
} else {
    optimiser <- match.arg(optimiser)
    hglmethod <- match.arg(hglmethod)
    rslt <- hmmUV(y,par0=par0,K=K,
                  method=method,hglmethod=hglmethod,optimiser=optimiser,
                  optimMethod=optimMethod,stationary=stationary,
                  mixture=mixture,cis=cis,tolerance=tolerance,
                  digits=digits,verbose=verbose,itmax=itmax,crit=crit,
                  bicm=bicm,X=X,addIntercept=addIntercept,
                  lmc=lmc,hessian=hessian,...)
}
# Tidy up the dimension names of tpm and Rho (or the
# levels of Rho$state, as is appropriate).
if(K > 1) {
    rownames(rslt$tpm) <- colnames(rslt$tpm) <- stnms
    if(cis) {
        names(rslt$ispd) <- stnms
    } else {
        if(length(y)==1) {
            rslt$ispd <- as.vector(rslt$ispd)
            names(rslt$ispd) <- stnms
        } else {
            colnames(rslt$ispd) <- stnms
        }
    }
    if(bivar) {
        if(indep) {
           colnames(rslt$Rho[[1]]) <- colnames(rslt$Rho[[2]]) <- stnms
        } else {
           dimnames(rslt$Rho)[[3]] <- stnms
        }
    } else {
        levels(rslt$Rho$state) <- stnms
        if(identical(prednames,"Intercept")) {
            rslt <- append(rslt,after=1,list(Rho.matrix=cnvrtRho(rslt$Rho)))
        }
    }
}

# Add the value of "stationary" to the result.
naft <- which(names(rslt)=="tpm")
rslt  <- append(rslt,list(stationary=stationary),after=naft)

# Add the lengths of the observation sequences to the returned value.
ylnths <- sapply(y,nrow)
if("prior.emsteps" %in% names(rslt)) {
    naft <- which(names(rslt)=="prior.emsteps")
} else {
    naft <- which(names(rslt)=="nstep")
}
rslt   <- append(rslt,list(ylengths=ylnths),after=naft)

# Add the missing value fraction(s) to the returned value.
naft   <- naft + 1
nafrac <- nafracCalc(y)
rslt   <- append(rslt,list(nafrac=nafrac),after=naft)
naft   <- naft + 1

# Possibly add y and X to the returned value.
if(keep.y) {
    rslt <- append(rslt,list(y=y),after=naft)
    naft <- naft + 1
}
keep.X <- keep.X & !bivar & !is.null(X)
if(keep.X) {
    rslt <- append(rslt,list(X=X),after=naft)
    naft <- naft + 1
}

# Add some auxiliary information to the returned value.
rslt <- append(rslt,list(parity=attr(y,"parity")),after=naft)
naft <- naft + 1
rslt <- append(rslt,list(numeric=attr(y,"numeric")),after=naft)
args <- list(method=method,optimiser=optimiser,
             optimMethod=optimMethod,stationary=stationary,
             mixture=mixture,cis=cis,tolerance=tolerance,
             itmax=itmax,crit=crit,addIntercept=addIntercept)
rslt$args   <- args
class(rslt) <- "hmm.discnp"
rslt
}
