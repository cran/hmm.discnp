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

# Check that one of par0 and K is specified.
if(is.null(par0) & is.null(K))
	stop('One of par0 and K must be specified.')

# Put together a data name tag for the output.
if(is.null(data.name)) data.name <- deparse(substitute(y))

# Transform the possible values of y to 1:nval where nval is
# the number of unique values of the original y.
y <- tidyup(y,yval)
nc <- ncol(y)
nval <- length(unique(y))

# If K=1 do the triv thing:
if(K==1) {
	Rho <- as.matrix(table(y)/length(y))
	ll  <- sum(log(ffun(y,Rho)))
	return(list(Rho=Rho,tpm=NA,ispd=NA,log.like=ll,
               converged=NA,nstep=NA,data.name=data.name))
}

# Pick out the index of the stopping criterion:
icrit <- match(crit,c('PCLL','L2','Linf'))
if(is.na(icrit)) stop('Stopping criterion not recognized.')

# Perform initial setting-up.
if(is.null(par0)) {
	if(is.null(rand.start)) rand.start <- list(tpm=FALSE,Rho=FALSE)
	par0  <- init.all(nval,K,rand.start,mixture)
}
else
	K     <- nrow(par0$tpm)

tpm  <- par0$tpm
ispd <- revise.ispd(tpm)
Rho  <- par0$Rho

m         <- nrow(Rho)
old.theta <- c(c(tpm[,-K]),c(Rho[1:(m-1),]))
old.ll    <- -Inf
digits    <- 2+ceiling(abs(log10(tolerance)))

if(verbose) cat('\n      Initial set-up completed ...\n\n')

# Set the level below which probabilities are considered to
# be noise:

# Update:
em.step <- 1
if(verbose) cat('Repeating ...\n\n')

chnge <- numeric(3)
repeat{
	if(verbose) cat(paste('EM step ',em.step,':\n',sep=''))

# Get the probabilities of the observations.
	fy <- ffun(y,Rho)

# Calculate the recursive probabilities.
	rp <- recurse(fy,tpm,ispd,nc)

# Calculate the parameters.
	tpm  <- revise.tpm(rp$xi,mixture)
	ispd <- revise.ispd(tpm)
	Rho  <- revise.rho(y,rp$gamma,m)

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
rp <- recurse(fy,tpm,ispd,nc)
ll <- sum(log(abs(rp$llc)))

# Return:
list(Rho=Rho,tpm=tpm,ispd=ispd,log.like=ll,converged=converged,
     nstep=nstep,data.name=data.name)
}
