hmm <- function(y,yval=NULL,par0=NULL,K=NULL,rand.start=NULL,
                mixture=FALSE,tolerance=1e-4,verbose=FALSE,itmax=200,
                crit='PCLL',data.name=NULL)
{
# Function hmm.  To fit a Hidden Markov model to a data set where the
# observations come from one of a number of finite discrete
# distributions, depending on the (hidden) state of the Markov chain.
# These distributions are specified by a matrix Rho = [rho_ij] where
# rho_ij = P(X = x_i | S = j), X being the observable random variable
# and S being the hidden state.
# 
# Note that y is allowed to be a matrix, each column being interpreted
# as an independent replicate of the observation sequence.
#
# Copyright (C) 1997 by T. Rolf Turner and Limin Liu.
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that both that copyright notice and this permission
# notice appear in supporting documentation.
#

# Check that one of par0 and K is specified.
if(is.null(par0) & is.null(K))
	stop('One of par0 and K must be specified.')

# Put together a data name tag for the output.
if(is.null(data.name)) data.name <- deparse(substitute(y))

# Pick out the index of the stopping criterion:
icrit <- match(crit,c('PCLL','L2','Linf'))
if(is.na(icrit)) stop('Stopping criterion not recognized.')

# Transform the possible values of y to 1:nval where nval is
# the number of unique values of the original y.
if(is.null(yval)) yval <- sort(unique(as.vector(y)))
nval <- length(yval)
if(is.matrix(y)) {
	nr <- nrow(y)
	nc <- ncol(y)
}
else {
	nr <- length(y)
	nc <- 1
}
y <- matrix(match(y,yval),nr,nc)

# Perform initial setting-up.
if(is.null(par0)) {
	if(is.null(rand.start)) rand.start <- list(tpm=FALSE,Rho=FALSE)
	par0  <- init.all(nval,K,rand.start,mixture)
}
else
	K     <- nrow(par0$tpm)

tpm <- par0$tpm
Rho <- par0$Rho

m         <- nrow(Rho)
old.theta <- c(c(tpm[,-K]),c(Rho[1:(m-1),]))
old.ll    <- -Inf
digits    <- 2+ceiling(abs(log10(tolerance)))

if(verbose) cat('\n      Initial set-up completed ...\n\n')

# Set the level below which probabilities are considered to
# be noise:
epsilon <- 10*.Machine$double.eps

# Update:
em.step <- 1
if(verbose) cat('Repeating ...\n\n')

chnge <- numeric(3)
repeat{
	if(verbose) cat(paste('EM step ',em.step,':\n',sep=''))

# Get the probabilities of the observations.
	fy <- ffun(y,Rho)

# Calculate the recursive probabilities.
	rp <- recurse(fy,tpm,nc,epsilon)

# Calculate the parameters.
	tpm <- revise.tpm(rp$xi,mixture)
	Rho <- revise.rho(y,rp$gamma,m)

# Test for convergence:
	new.theta <- c(c(tpm[,-K]),c(Rho[1:(m-1),]))
	ll <-  sum(log(abs(rp$llc)))
	chnge[1] <- if(old.ll > -Inf)
			100*(ll - old.ll)/abs(old.ll)
		    else
			Inf
	chnge[2] <- sqrt(sum((old.theta-new.theta)**2))
	chnge[3] <- max(abs(old.theta-new.theta))
	if(verbose){
		cat('     Log-likelihood: ',
            	format(round(ll,digits)),'\n',sep='')
		cat('     Percent decrease in log-likelihood: ',
		    format(round(chnge[1],digits)),'\n',sep='')
		cat('     Root-SS of change in coef.: ',
		    format(round(chnge[2],digits)),'\n',sep='')
		cat('     Max. abs. change in coef.: ',
		    format(round(chnge[3],digits)),'\n',sep='')
	}

	if(chnge[icrit] < tolerance) {
			theta <- new.theta
			converged <- TRUE
			nstep <- em.step
			break
		}

	if(em.step >= itmax) {
		cat('Failed to converge in ',itmax,' EM steps.\n',sep='')
		theta <- new.theta
		converged <- FALSE
		nstep <- em.step
		break
	}
	old.theta <- new.theta
	old.ll    <- ll
	em.step   <- em.step + 1
}

# Calculate the log-likelihood of the current (final)
# parameter values.
fy <- ffun(y,Rho)
rp <- recurse(fy,tpm,nc,epsilon)
ll <- sum(log(abs(rp$llc)))

# Calculate the initial state probability distribution.
ispd   <- revise.ispd(tpm)

# Return:
list(Rho=Rho,tpm=tpm,ispd=ispd,log.like=ll,converged=converged,
     nstep=nstep,data.name=data.name)

}
