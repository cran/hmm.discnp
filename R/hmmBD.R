hmmBD <- function(y,par0,K,stationary,
                  mixture,cis,tolerance,digits,
                  verbose,itmax,crit,bicm) {
#
# Function hmmBD.  To conduct the fitting of a Hidden Markov model
# when the observations are bivariate and NOT assumed to be conditionally
# independent.

# Check for adequacy of data.
yawl <- do.call(rbind,y)
X    <- yawl[,1]
Y    <- yawl[,2]
Tab  <- table(X,Y)
if(sum(Tab)==0) {
    whinge <- paste0("No non-missing correspondences between variable 1 and\n",
                     "variable 2.  The data are inadequate for fitting a\n",
                     "bivariate dependent model.\n")
    stop(whinge)
}

# If K=1 do the triv thing.  (Not quite so triv in the
# bivariate dependent setting! Handling missing values is
# much more complicated than in the univariate or bivariate
# independent settings.)
if(K==1) {
    lvls <- attr(y,"lvls")
    dnms <- c(lvls,list("1"))
    ym   <- do.call(rbind,y)
    X    <- factor(ym[,1],levels=c(lvls[[1]],NA),exclude=NULL)
    Y    <- factor(ym[,2],levels=c(lvls[[2]],NA),exclude=NULL)
    G    <- table(X,Y,useNA="always")
    m    <- nrow(G)
    n    <- ncol(G)
    Rho0 <- G[-m,-n]
    Rho0 <- Rho0/sum(Rho0)
    Rho0 <- array(Rho0,dim=c(dim(Rho0),1))
    dimnames(Rho0) <- dnms
    G    <- array(G,dim=c(dim(G),1))
    Rho  <- msRho(Rho0,G)
    dimnames(Rho) <- dnms
    ll   <- sum(log(ffun(y,Rho,type=3)))
    npar <- prod(dim(Rho))-1
    AIC  <- -2*ll+2*npar
    BIC  <- -2*ll+bicm*npar
    rslt <- list(Rho=Rho,tpm=NA,ispd=NA,log.like=ll,par0=NA,npar=npar,
                 converged=NA,nstep=NA,
                 stationary=NA,cis=NA,AIC=AIC,BIC=BIC)
    class(rslt) <- "hmm.discnp"
    return(rslt)
}

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
Rho    <- par0$Rho
m1     <- dim(Rho)[1]
m2     <- dim(Rho)[2]

# Get the lengths of the observations.
lns  <- sapply(y,nrow)

# Set the number of digits with which to print out
# "progress reports".
if(is.null(digits)) digits <- 2+ceiling(abs(log10(tolerance)))

eyedrop   <- cumsum(rep(prod(dim(Rho)[-3]),dim(Rho)[3]))
old.theta <- c(as.vector(tpm[,-K]),as.vector(Rho)[-eyedrop])
fy        <- ffun(y,Rho,type=3)
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
	Rho  <- revise.rho(y,rp$gamma,type=3)

# Update the log likelihood on the basis of the
# new parameter estimates.  This entails calculating
# the new recursive probabilities (which will be used
# to update the parameter estimates on the *next* EM
# step, if necessary).
	fy <- ffun(y,Rho,type=3)
	rp <- recurse(fy,tpm,ispd,lns)
	ll <-  sum(log(rp$llc))

# Test for convergence:
        new.theta <- c(as.vector(tpm[,-K]),as.vector(Rho)[-eyedrop])
	chnge[1]  <- 100*(ll - old.ll)/(abs(old.ll) + tolerance)
	chnge[2]  <- sqrt(sum((old.theta-new.theta)**2))/
                     (sqrt(sum(new.theta)^2) + tolerance)
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
npar <- nispar + K*(K-1) + prod(dim(Rho))-K

AIC  <- -2*ll + 2*npar
BIC  <- -2*ll+bicm*npar
names(chnge) <- c("PCCL","L1","Linf")
rslt <- list(Rho=Rho,tpm=tpm,ispd=ispd,log.like=ll,stopCrit=chnge,par0=par0,
             npar=npar,bicm=bicm,stopCrit=chnge,converged=converged,nstep=nstep,
             AIC=AIC,BIC=BIC)
rslt
}
