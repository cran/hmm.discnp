hmmBI <- function(y,yval=NULL,par0=NULL,K=NULL,rand.start=NULL,stationary=cis,
                  mixture=FALSE,cis=TRUE,tolerance=NULL,digits=NULL,verbose=FALSE,
                  itmax=200,crit,bicm) {
#
# Function hmmBI.  To conduct the fitting of a Hidden Markov model
# when the observations are bivariate and assumed to be conditionally
# independent.

# Check that the observation values are compatible
# with yval if it is specified.

# If Rho has been specified in par0, and if its components have
# row names, replace the appropriate components of yval by these
# row names.
if(!is.null(par0$Rho)) {
    satch <- vector("list",2)
    rval  <- rep(FALSE,2)
    for(i in 1:2) {
        rnms <- rownames(par0$Rho[[i]])
        if(!is.null(rnms)) {
            satch[[i]] <- rnms
            rval[i] <- TRUE
        }
    }
    if(all(rval)) {
        if(!is.null(yval)) {
            whinge <- paste0("The specified initial value of Rho has components\n",
                             "with row names.  These take precedence and yval has\n",
                             "been replaced by these row names.\n")
            warning(whinge)
        } else if(sum(rval) == 1) {
            whinge <- paste0("Inconsistency in the presence of row names in the\n",
                             "components of \"Rho\".\n")
            stop(whinge)
        }
        yval <- satch
    }
}

lvls <- attr(y,"lvls")
if(is.null(yval)) yval <- lvls
if(!(all(lvls[[1]]%in%yval[[1]]) & all(lvls[[2]]%in%yval[[2]])))
        stop("Specified y values do not include all observed y values.\n")
nval <- sapply(yval,length)
lns  <- sapply(y,nrow)

# Now replace the "lvls" attribute of "y" by "yval" so that this
# attribute really does provide the appropriate levels.
attr(y,"lvls") <- yval

# If K=1 do the triv thing:
if(K==1) {
        Dat <- makeDat(y,X=NULL)
        Rho <- vector("list",2)
        DF  <- do.call(rbind,Dat)
	Rho[[1]] <- as.matrix(table(DF[,1]))
        Rho[[1]] <- Rho[[1]]/sum(Rho[[1]])
        rownames(Rho[[1]]) <- yval[[1]]
	Rho[[2]] <- as.matrix(table(DF[,2]))
        Rho[[2]] <- Rho[[2]]/sum(Rho[[2]])
        rownames(Rho[[2]]) <- yval[[2]]
	ll   <- sum(log(ffun(Dat,Rho,type=3)))
        npar <- sum(sapply(Rho,nrow))-2
        AIC  <- -2*ll+2*npar
        BIC  <- -2*ll+bicm*npar
        rslt <- list(Rho=Rho,tpm=NA,ispd=NA,log.like=ll,par0=NA,npar=npar,
                     converged=NA,nstep=NA,
                     stationary=NA,cis=NA,AIC=AIC,BIC=BIC)
        class(rslt) <- "hmm.discnp"
        return(rslt)
}

# If par0 was not specified, initialise it.
if(is.null(par0)) {
	if(is.null(rand.start)) rand.start <- list(tpm=FALSE,Rho=FALSE)
	par0 <- init.all(nval,K,rand.start,mixture,indep=TRUE)
} else {
    K <- nrow(par0$tpm)
    if(K != ncol(par0$tpm))
            stop("The specified tpm is not square.\n")
}

for(i in 1:2) {
    ordnl <- if(i==1) "first" else "second"
    if(nrow(par0$Rho[[i]]) < nval[i])
        stop(paste("The row dimension of Rho[[",i,"]] is less than\n",
                   "the number of distinct values of the ",ordnl," variate.\n"))
}
rnms <- vector("list",2)
for(i in 1:2) {
    if(is.null(rownames(par0$Rho[[i]]))) {
        if(nval[i] != nrow(par0$Rho[[1]])) {
            whinge <- paste0("No rownames for Rho[[",i,"]] and nrow(Rho[[",i,"]])",
                             " is not\n",
                             "equal to the number of unique values of the first",
                             " variate.\n")
            stop(whinge)
        }
        rnms[[i]] <- rownames(par0$Rho[[i]]) <- yval[[i]]
    } else {
        rnms[[i]] <- rownames(par0$Rho[[i]])
        if(!all(yval[[i]]%in%rnms[[i]])) {
            whinge <- paste0("The row names of the initial value of Rho[[",i,"]]",
                             " do not\n",
                            "include all possible values of the first variate.\n")
            stop(whinge)
        }
    }
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
Rho    <- par0$Rho
m1     <- nrow(Rho[[1]])
m2     <- nrow(Rho[[2]])

# Set the tolerance and the number of digits with which to print out
# "progress reports".
if(is.null(tolerance)) tolerance <- 1e-4
if(is.null(digits)) digits <- 2+ceiling(abs(log10(tolerance)))

old.theta <- c(as.vector(tpm[,-K]),as.vector(Rho[[1]][1:(m1-1),]),
                                   as.vector(Rho[[2]][1:(m2-1),]))
fy        <- ffun(Dat,Rho,type=3)
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
	fy <- ffun(Dat,Rho,type=3)
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

# Tidy up a bit:
if(length(y)==1) ispd <- as.vector(ispd)
stnms <- rownames(par0$tpm)
if(is.null(stnms)) stnms <- 1:K
rownames(tpm) <- stnms
colnames(tpm) <- stnms
names(ispd)   <- stnms
colnames(Rho[[1]]) <- stnms
colnames(Rho[[2]]) <- stnms

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
