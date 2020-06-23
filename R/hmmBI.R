hmmBI <- function(y,par0,K,stationary,
                  mixture,cis,tolerance,digits,verbose,
                  itmax,crit,bicm) {
#
# Function hmmBI.  To conduct the fitting of a Hidden Markov model
# when the observations are bivariate and assumed to be conditionally
# independent.  Called by hmm().  Not intended to be called directly.

# If K=1 do the triv thing:
if(K==1) {
        lvls <- attr(y,"lvls")
        Dat <- makeDat(y,X=NULL)
        Rho <- vector("list",2)
        DF  <- do.call(rbind,Dat)
	Rho[[1]] <- as.matrix(table(DF[,1]))
        Rho[[1]] <- Rho[[1]]/sum(Rho[[1]])
        rownames(Rho[[1]]) <- lvls[[1]]
	Rho[[2]] <- as.matrix(table(DF[,2]))
        Rho[[2]] <- Rho[[2]]/sum(Rho[[2]])
        rownames(Rho[[2]]) <- lvls[[2]]
	ll   <- sum(log(ffun(Dat,Rho,type=2)))
        npar <- sum(sapply(Rho,nrow))-2
        AIC  <- -2*ll+2*npar
        BIC  <- -2*ll+bicm*npar
        rslt <- list(Rho=Rho,tpm=NA,ispd=NA,log.like=ll,par0=NA,npar=npar,
                     converged=NA,nstep=NA,
                     stationary=NA,cis=NA,AIC=AIC,BIC=BIC)
        class(rslt) <- "hmm.discnp"
        return(rslt)
}

# Turn the observation data into a list of *data frames* whose
# columns are *factors*.
Dat  <- makeDat(y,X=NULL)

# Pick out the index of the stopping criterion:
icrit <- match(crit,c('PCLL','L2','Linf'))
if(is.na(icrit)) stop(paste("Stopping criterion",crit,"not recognized.\n"))

# Perform initial setting-up.
tpm    <- par0$tpm
if(stationary) {
    ispd   <- revise.ispd(tpm)
} else { # Make the chains equally likely to start in any state.
    ispd <- matrix(1/K,K,length(y))
}
Rho   <- par0$Rho
m1    <- nrow(Rho[[1]])
m2    <- nrow(Rho[[2]])
lns   <- sapply(y,nrow)

# Set the tolerance and the number of digits with which to print out
# "progress reports".
if(is.null(tolerance)) tolerance <- 1e-4
if(is.null(digits)) digits <- 2+ceiling(abs(log10(tolerance)))

old.theta <- c(as.vector(tpm[,-K]),as.vector(Rho[[1]][1:(m1-1),]),
                                   as.vector(Rho[[2]][1:(m2-1),]))
fy        <- ffun(Dat,Rho,type=2)
rp        <- recurse(fy,tpm,ispd,lns)
old.ll    <- sum(log(rp$llc))

# Ready to go.
    if(verbose){
        cat("\n      Initial set-up completed ...\n")
        cat("\n      Initial log-likelihood: ",
        format(round(old.ll,digits)),"\n\n",sep="")
    }

# Revise:
em.step <- 1
if(verbose) cat('Repeating ...\n\n')
chnge <- numeric(3)
repeat{
	if(verbose) cat(paste('EM step ',em.step,':\n',sep=''))

# Calculate the parameters.
	tpm  <- revise.tpm(rp$xi,mixture)
	ispd <- if(stationary) {
			revise.ispd(tpm)
		} else {
			revise.ispd(gamma=rp$gamma,lns=lns,cis=cis)
		}
	Rho  <- revise.rho(y,rp$gamma,type=2)

# Update the log likelihood on the basis of the
# new parameter estimates.  This entails calculating
# the new recursive probabilities (which will be used
# to update the parameter estimates on the *next* EM
# step, if necessary).
	fy <- ffun(Dat,Rho,type=2)
	rp <- recurse(fy,tpm,ispd,lns)
	ll <-  sum(log(rp$llc))

# Test for convergence:
	new.theta <- c(as.vector(tpm[,-K]),as.vector(Rho[[1]][1:(m1-1),]),
                                           as.vector(Rho[[2]][1:(m2-1),]))
	chnge[1]  <- 100*(ll - old.ll)/(abs(old.ll) + tolerance)
	chnge[2]  <- sqrt(sum((old.theta-new.theta)^2))/
                     (sqrt(sum(new.theta^2)) + tolerance)
	chnge[3]  <- max(abs(old.theta-new.theta))/
                     (max(abs(new.theta)) + tolerance)
	if(verbose){
		cat('     Log-likelihood: ',
            	format(round(ll,digits)),'\n',sep='')
		cat('     Scaled percent increase in log-likelihood: ',
		    format(round(chnge[1],digits)),'\n',sep='')
		cat('     Scaled root-SS of change in coef.: ',
		    format(round(chnge[2],digits)),'\n',sep='')
		cat('     Scaled max. abs. change in coef.: ',
		    format(round(chnge[3],digits)),'\n',sep='')
	}

	if(chnge[icrit] < tolerance) {
			converged <- TRUE
			nstep <- em.step
			break
		}

	if(em.step >= itmax) {
		cat('Failed to converge in ',itmax,' EM steps.\n',sep='')
		converged <- FALSE
		nstep <- em.step
		break
	}

# Replace the ``old'' parameter and log likelihood values
# by the new ones.
	old.theta <- new.theta
	old.ll    <- ll
# Increment the step number.
	em.step   <- em.step + 1
}

if(stationary) {
    nispar <- 0
} else {
    nispar <- if(cis) K-1 else (K-1)*length(lns)
}
npar <- nispar + K*(K-1) + K*(sum(sapply(Rho,nrow))-2)

AIC  <- -2*ll + 2*npar
BIC  <- -2*ll+bicm*npar
names(chnge) <- c("PCCL","L1","Linf")

list(Rho=Rho,tpm=tpm,ispd=ispd,log.like=ll,stopCrit=chnge,par0=par0,
             npar=npar,bicm=bicm,converged=converged,nstep=nstep,
             AIC=AIC,BIC=BIC)
}
