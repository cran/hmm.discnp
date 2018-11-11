hmmUV <- function(y,yval=NULL,par0=NULL,K=NULL,rand.start=NULL,
                  method=c("EM","bf","LM","SD"),hglmethod=c("fortran","oraw","raw"),
                  optimiser=c("nlm","optim"),optimMethod=NULL,stationary=cis,
                  mixture=FALSE,cis=TRUE,tolerance,digits=NULL,verbose=FALSE,
                  itmax=200,crit,bicm,newstyle=TRUE,X=NULL,addIntercept=TRUE,
                  lmc=10,hessian=FALSE,...) {
#
# Function hmmUV.  To conduct the fitting of a Hidden Markov model
# when the observations are univariate.  Note that in the univariate
# case we are allowed to have a matrix X of auxiliary predictors
# if the method is "EM".  Or "bf"???
#

# Check that the newstyle value is legitimate.
if(!(is.null(X) | newstyle))
    stop("If \"X\" is non-NULL then \"newstyle\" must be TRUE.\n")

# Check that the observation values are compatible
# with yval if it is specified.
lvls <- attr(y,"lvls")
if(is.null(yval)) yval <- lvls
if(!all(lvls%in%yval))
        stop("Specified y values do not include all observed y values.\n")
lns  <- sapply(y,length)
nval <- length(yval)

# Now replace the "lvls" attribute of "y" by "yval" so that this
# attribute really does provide the appropriate levels.
attr(y,"lvls") <- yval

# Make sure that X, if supplied, has the appropriate form.
if(!is.null(X)) {
    X <- tidyList(X,rp="predictor",addIntercept=addIntercept)
    checkyXoK(y,X)
    prednames <- attr(X,"prednames")
} else prednames <- "Intercept"

# If K=1 do the triv thing.  (Not so triv anymore, what with "newstyle"
# and "X" rearing their pretty little heads!) If X is provided, call
# upon multinom to do the estimation of Rho.

if(K==1) {
    yf <- factor(unlist(y),levels=yval)
    if(newstyle) {
        if(is.null(X)) {
           XX <- as.data.frame(matrix(1,nrow=length(yf),ncol=1))
           pnms <- "Intercept"
        } else {
           XX <- as.data.frame(do.call(rbind,X))
           pnms <- attr(X,"prednames")
        }
        mnfit <- multinom(yf ~ 0 + ., data=XX,trace=FALSE)
        CCC   <- rbind(coef(mnfit),0)
        colnames(CCC) <- pnms
        rownames(CCC) <- 1:nrow(CCC)
        Rho   <- data.frame(y=factor(yval),state=factor(rep(1,length(yval))),CCC)
        npar  <- (nrow(Rho) - 1)*(ncol(Rho)-2)
    } else {
        ttt   <- table(yf)
        Rho   <- as.matrix(ttt/sum(ttt))
        npar  <- nrow(Rho) - 1
    }
    Dat  <- makeDat(y,X,addIntercept=addIntercept)
    ll   <- sum(log(ffun(Dat,Rho,type=if(newstyle) 1 else 2)))
    AIC  <- -2*ll+2*npar
    BIC  <- -2*ll+bicm*npar
    rslt <- list(Rho=Rho,tpm=NA,ispd=NA,log.like=ll,par0=NA,
                 npar=npar,converged=NA,nstep=NA,AIC=AIC,BIC=BIC)
    class(rslt) <- "hmm.discnp"
    return(rslt)
}

# If par0 was not specified, initialise it.
if(is.null(par0)) {
    if(is.null(rand.start)) rand.start <- list(tpm=FALSE,Rho=FALSE)
    par0  <- init.all(nval,K,rand.start,mixture,indep=NA,newstyle,yval,
                      prednames)
}
else {
    if(newstyle) {
        if(!is.data.frame(par0$Rho))
            stop("Component \"Rho\" of \"par0\" is of the wrong form.\n")
        if(length(levels(par0$Rho$y)) < nval)
            stop(paste("The number of levels of par0$Rho$y is less than\n",
                       "the number of distinct y-values.\n"))
        if(length(levels(par0$Rho$state)) != K)
            stop("The number of levels of par0$Rho$state is not equal to K.\n")
    } else {
        if(nrow(par0$Rho) < nval)
            stop(paste("The row dimension of Rho is less than\n",
                       "the number of distinct y-values.\n"))
    }
}

if(newstyle) {
    rhovals <- levels(par0$Rho$y)
    if(!all(yval%in%rhovals)) {
        whinge <- paste("The levels of Rho$y in the initial value of \"Rho\"\n",
                        "do not include all possible y-values.\n")
        stop(whinge)
    }
} else {
    if(is.null(rownames(par0$Rho))) {
        if(length(yval) != nrow(par0$Rho)) {
            whinge <- paste("No rownames for Rho and nrow(Rho) is not equal\n",
                            "to the number of unique y values.\n",sep="")
            stop(whinge)
        }
        rhovals <- rownames(par0$Rho) <- yval
    } else {
        rhovals <- rownames(par0$Rho)
        if(!all(yval%in%rhovals)) {
            whinge <- paste("The row names of the initial value of \"Rho\" do not\n",
                        "include all possible y-values.\n")
            stop(whinge)
        }
    }
}

# Add ispd to par0 if it isn't there already.
ispd <- par0$ispd
if(is.null(ispd)) {
    if(stationary) {
        ispd <- revise.ispd(par0$tpm)
    } else { # Make the chains equally likely to start in any state.
        ispd <- if(cis) rep(1/K,K) else matrix(1/K,K,length(y))
    }
    par0 <- c(list(ispd=ispd),par0)
}

# Set the parameter count.
if(stationary) {
    nispar <- 0
} else {
    nispar <- if(cis) K-1 else (K-1)*length(lns)
}
nrhopar <- if(newstyle)
    (nrow(par0$Rho) - K)*(ncol(par0$Rho)-2) else K*(nrow(par0$Rho) - 1)
ntpmpar <- K*(K-1)
npar    <- nispar + ntpmpar + nrhopar

# If newstyle is TRUE and method is either EM or bf, put together
# a list of data frames consisting of y and the predictors.
# Otherwise make Dat a list of factors.

Dat  <- makeDat(y,X,addIntercept=addIntercept)

# Choose method.
if(method=="EM") {

# Pick out the index of the stopping criterion:
    icrit <- match(crit,c("PCLL","L2","Linf"))

# Perform initial setting-up.
    tpm       <- par0$tpm
    ispd      <- par0$ispd
    Rho       <- par0$Rho
    t1        <- as.vector(tpm[,-K])
    t2        <- paramExtract(Rho,newstyle)
    old.theta <- c(t1,t2)
    fy        <- ffun(Dat,Rho,type=if(newstyle) 1 else 2)
    rp        <- recurse(fy,tpm,ispd,lns)
    old.ll    <- sum(log(rp$llc))

# Set the number of digits with which to print out
# "progress reports".
if(is.null(digits)) digits <- 2+ceiling(abs(log10(tolerance)))

# Ready to go.
    if(verbose){
        cat("\n      Initial set-up completed ...\n")
        cat("\n      Initial log-likelihood: ",
        format(round(old.ll,digits)),"\n\n",sep="")
    }
    
# Revise:
    em.step <- 1
    if(verbose) cat("Repeating ...\n\n")
    chnge <- numeric(3)
    repeat{
    	if(verbose) cat(paste("EM step ",em.step,":\n",sep=""))
    
# Calculate the parameters.
    	tpm  <- revise.tpm(rp$xi,mixture)
    	ispd <- if(stationary) {
    			revise.ispd(tpm)
    		} else {
    			revise.ispd(gamma=rp$gamma,lns=lns,cis=cis)
    		}
    	Rho  <- revise.rho(Dat,rp$gamma,type=if(newstyle) 1 else 2)
    
# Update the log likelihood on the basis of the
# new parameter estimates.  This entails calculating
# the new recursive probabilities (which will be used
# to update the parameter estimates on the *next* EM
# step, if necessary).
    	fy <- ffun(Dat,Rho,type=if(newstyle) 1 else 2)
    	rp <- recurse(fy,tpm,ispd,lns)
    	ll <-  sum(log(rp$llc))

# Previously the following line was
#       chnge[1]  <- 100*(ll - old.ll)/abs(old.ll)
# Changed it (30/06/2018) for consistency with hmmLM().
    	chnge[1]  <- 100*(ll - old.ll)/(abs(old.ll + tolerance))
        if(ll < old.ll) {
            Xused <- !is.null(X)
            if(verbose) {
                cat("Decrease in log likelihood --- by ",abs(chnge[1]),
                    " percent.\n",sep="")
                cat("Switching to the", if(Xused) "\"bf\"" else "\"LM\"","method.\n")
            }
            par  <- list(ispd=ispd,tpm=tpm,Rho=Rho)
            if(Xused) {
                rslt <- hmmNumOpt(Dat,par0=par,stationary=stationary,verbose=verbose,
                                  itmax=itmax,bicm=bicm,rhovals=rhovals,npar=npar,
                                  optimiser=optimiser,optimMethod=optimMethod,
                                  hessian=hessian,...)
                method <- "bf"
            } else {
                Dat  <- lapply(Dat,function(x){x$y})
                rslt <- hmmLM(Dat,par0=par,itmax=itmax,crit=crit,lmc=lmc,
                              tolerance=tolerance,bicm=bicm,rhovals=rhovals,
                              hglmethod=hglmethod,digits=digits,verbose=verbose)
                method <- "LM"
            }
            rslt$prior.emsteps <- em.step
            break
        } 
# Test for convergence:
        t1        <- as.vector(tpm[,-K])
        t2        <- paramExtract(Rho,newstyle)
        new.theta <- c(t1,t2)
    	chnge[2]  <- sqrt(sum((old.theta-new.theta)**2))/
                     (sqrt(sum(new.theta^2)) + tolerance)
    	chnge[3]  <- max(abs(old.theta-new.theta))/
                     (max(abs(new.theta)) + tolerance)
    	if(verbose){
    	    cat("     Log-likelihood: ",
                format(round(ll,digits)),"\n",sep="")
    	    cat("     Scaled percent increase in log-likelihood: ",
                format(round(chnge[1],digits)),"\n",sep="")
    	    cat("     Scaled root-SS of change in coef.: ",
                format(round(chnge[2],digits)),"\n",sep="")
    	    cat("     Scaled max. abs. change in coef.: ",
                format(round(chnge[3],digits)),"\n",sep="")
    	}
    
    	if(chnge[icrit] < tolerance) {
    			converged <- TRUE
    			nstep <- em.step
    			break
    		}
    
    	if(em.step >= itmax) {
    		if(verbose) cat("Failed to converge in ",itmax," EM steps.\n",sep="")
    		converged <- FALSE
    		nstep <- em.step
    		break
    	}
    
# Replace the "old" parameter and log likelihood values
# by the new ones.
    	old.theta <- new.theta
    	old.ll    <- ll
# Increment the step number.
    	em.step   <- em.step + 1
    }

if(method=="EM") { # Might have switched to "bf" or "LM" due
                   # to decrease in log likelihood.
# Tidy up ispd.
    if(length(y)==1) ispd <- as.vector(ispd)

# Set AIC and BIC.
    AIC     <- -2*ll + 2*npar
    BIC     <- -2*ll + bicm*npar

# Put together the result.
    names(chnge) <- c("PCCL","L1","Linf")
    rslt    <- list(Rho=Rho,tpm=tpm,ispd=ispd,log.like=ll,stopCrit=chnge,
                    par0=par0,npar=npar,bicm=bicm,converged=converged,nstep=nstep,
                    AIC=AIC,BIC=BIC)
}

} else {
    if(!stationary) {
        whinge <- paste0("None of the Levinberg-Marquardt algorithm, brute\n",
                         "force maximisation, or steepest descent methods",
                         " can\n", "(currently) handle non-stationary models.\n")
        stop(whinge)
    }
    if(method=="bf") {
        rslt <- hmmNumOpt(Dat,par0=par0,stationary=stationary,verbose=verbose,
                          itmax=itmax,bicm=bicm,rhovals=rhovals,npar=npar,
                          optimiser=optimiser,optimMethod=optimMethod,
                          hessian=hessian,...)
        rslt$prior.emsteps <- 0
    } else if(method=="LM") {
        if(!is.null(X)) {
            whinge <- paste0("Cannot currently use the Levenberg-Marquardt",
                             " algorithm when\n there are auxiliary",
                             " predictors (\"X\").\n") 
            stop(whinge)
        }
        rslt <- hmmLM(Dat,par0=par0,itmax=itmax,crit=crit,lmc=lmc,bicm=bicm,
                      rhovals=rhovals,hglmethod=hglmethod,
                      tolerance=tolerance,digits=digits,verbose=verbose)
        rslt$prior.emsteps <- 0
    } else if(method=="SD") {
        if(!is.null(X)) {
            whinge <- paste0("Cannot currently use the method of steepest",
                             " descent when\n there are auxiliary",
                             " predictors (\"X\").\n") 
            stop(whinge)
        }
        rslt <- hmmSD(Dat,par0=par0,itmax=itmax,crit=crit,bicm=bicm,
                      rhovals=rhovals,hglmethod=hglmethod,tolerance=tolerance,
                      digits=digits,verbose=verbose)
    } else stop(paste0("Method ",method," not recognised.\n"))
}

# Assign column names (names of the states) to the tpm and Rho
# components of "rslt" and names to the ispd component.
stnms <- rownames(par0$tpm)
if(!is.null(stnms)) {
    rownames(rslt$tpm) <- stnms
    colnames(rslt$tpm) <- stnms
    if(is.matrix(rslt$ispd)) {
        colnames(rslt$ispd) <- stnms
    } else {
        names(rslt$ispd)   <- stnms
    }
    if(!newstyle) colnames(rslt$Rho) <- stnms
}
return(rslt)
}
