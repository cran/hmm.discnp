hmmUV <- function(y,par0,K,method,
                  hglmethod,optimiser,optimMethod,stationary,
                  mixture,cis,tolerance,digits,verbose,
                  itmax,crit,bicm,X,addIntercept,lmc,hessian,...) {
#
# Function hmmUV.  To conduct the fitting of a Hidden Markov model
# when the observations are univariate.  Note that in the univariate
# case we are allowed to have a matrix X of auxiliary predictors
# if the method is "EM" or "bf" (or if K=1).
#

# If X was supplied, get the predictor names.
if(!is.null(X)) {
    prednames <- attr(X,"prednames")
} else prednames <- "Intercept"

# If K=1 do the triv thing.  (Not so triv anymore, what with the
# data frame representation of Rho, and possible predictor variables
# X rearing their pretty little heads!) If X is provided, call
# upon multinom to do the estimation of Rho.

if(K==1) {
    lvls <- attr(y,"lvls")
    yf <- factor(unlist(y),levels=lvls)
    if(is.null(X)) {
        ttt    <- table(yf)
        caviar <- as.matrix(ttt/sum(ttt))
        Rho    <- cnvrtRho(caviar)
    } else {
        XX <- as.data.frame(do.call(rbind,X))
        mnfit <- multinom(yf ~ 0 + ., data=XX,trace=FALSE)
        CCC   <- rbind(coef(mnfit),0)
        colnames(CCC) <- prednames
        rownames(CCC) <- 1:nrow(CCC)
        Rho  <- data.frame(y=factor(lvls),state=factor(rep(1,length(lvls))),CCC)
    }
    npar <- (nrow(Rho) - 1)*(ncol(Rho)-2)
    Dat  <- makeDat(y,X,addIntercept=addIntercept)
    ll   <- sum(log(ffun(Dat,Rho,type=1)))
    AIC  <- -2*ll+2*npar
    BIC  <- -2*ll+bicm*npar
    rslt <- list(Rho=Rho,tpm=NA,ispd=NA,log.like=ll,par0=NA,
                 npar=npar,converged=NA,nstep=NA,AIC=AIC,BIC=BIC)
    class(rslt) <- "hmm.discnp"
    return(rslt)
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
lns <- sapply(y,length)
if(stationary) {
    nispar <- 0
} else {
    nispar <- if(cis) K-1 else (K-1)*length(lns)
}
nrhopar <- (nrow(par0$Rho) - K)*(ncol(par0$Rho)-2)
ntpmpar <- K*(K-1)
npar    <- nispar + ntpmpar + nrhopar
Rho     <- par0$Rho
rhovals <- levels(Rho$y)

# Put together a list of data frames consisting of y and the
# predictors.
Dat  <- makeDat(y,X,addIntercept=addIntercept)

# Choose method.
if(method=="EM") {

# Pick out the index of the stopping criterion:
    icrit <- match(crit,c("PCLL","L2","Linf"))

# Perform initial setting-up.
    tpm       <- par0$tpm
    ispd      <- par0$ispd
    t1        <- as.vector(tpm[,-K])
    t2        <- paramExtract(Rho)
    old.theta <- c(t1,t2)
    fy        <- ffun(Dat,Rho,type=1)
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
        Rho  <- revise.rho(Dat,gamma=rp$gamma,type=1)

# Update the log likelihood on the basis of the
# new parameter estimates.  This entails calculating
# the new recursive probabilities (which will be used
# to update the parameter estimates on the *next* EM
# step, if necessary).
        fy <- ffun(Dat,Rho,type=1)
    	rp <- recurse(fy,tpm,ispd,lns)
    	ll <-  sum(log(rp$llc))

    	chnge[1]  <- 100*(ll - old.ll)/(abs(old.ll + tolerance))
        if(ll < old.ll) {
            decll  <- round(abs(chnge[1]),digits)
            whinge <- paste0("\n  At the ",em.step,"th step, there was a DECREASE\n",
                             "  in the log likelihood in the amount of ",decll,".\n",
                             "  The EM algorithm appears not to be working.\n",
                             "  See the help for an explanation.  Change to one\n",
                             "  of the other methods, perhaps using the current\n",
                             "  results as starting values.\n")

            message(whinge)
            converged <- FALSE
            nstep <- em.step-1
            break
        }
# Test for convergence:
        t1        <- as.vector(tpm[,-K])
        t2        <- paramExtract(Rho)
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

rslt
}
