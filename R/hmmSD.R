hmmSD <- function(y,par0,itmax=200,crit,tolerance,
                  bicm,rhovals,hglmethod,digits=NULL,verbose=FALSE) {
#
# Function hmmSD.  The engine to carry out the fitting of a
# hidden Markov model with discrete emissions, modelled
# non-parametrically, using the method of steepest descent.
#

# Do some initial housekeeping.
K      <- nrow(par0$tpm)
theta  <- reparam(par0,stationary=TRUE)
npar   <- length(theta)
npro   <- K*(K-1)
old.ll <- get.l(theta,K,y)
oath   <- theta
if(is.null(tolerance)) tolerance <- 1e-4
if(is.null(digits)) digits <- 2+ceiling(abs(log10(tolerance)))
scrit  <- numeric(4)
icrit  <- match(crit,c("PCLL","L2","Linf","ABSGRD"))
if(is.na(icrit)) stop(paste("Stopping criterion",crit,"not recognized.\n"))

# Ready to go.
    if(verbose){
        cat("\n      Initial set-up completed ...\n")
        cat("\n      Initial log-likelihood: ",
        format(round(old.ll,digits)),"\n\n",sep="")
    }

# Update:
if(verbose) cat('Repeating ...\n\n')
sd.step <- 1

repeat { # Swear again; i.e. recurse.
	if(verbose) cat(paste('SD step ',sd.step,':\n',sep=''))
	theta <- steepest(K,y,theta)
        tpm   <- getTpm(theta,K,TRUE)
	if(identical(all.equal(tpm,diag(K)),TRUE)) {
		return(list(converged=FALSE,message="tpm equals identity"))
	}
        xxx   <- get.gl(theta,K,y)
	ll    <- xxx$ll
	grad  <- xxx$grad
        # The "old.ll > -Inf" test should no longer be necessary
        # but it does no real harm.
	scrit[1] <- if(old.ll > -Inf)
                        100*(ll - old.ll)/(abs(ll) + tolerance)
                    else Inf
        scrit[2] <- sqrt(sum((theta-oath)^2))/(sqrt(sum(theta^2)) + tolerance)
        scrit[3] <- max(abs(theta-oath))/(max(abs(theta)) + tolerance)
        scrit[4] <- sum(abs(grad))
        if(verbose){
                cat('     Log-likelihood: ',
                format(round(ll,digits)),'\n',sep='')
                cat('     Scaled percent increase in log-likelihood: ',
                    format(round(scrit[1],digits)),'\n',sep='')
                cat('     Scaled root-SS of change in coef.: ',
                    format(round(scrit[2],digits)),'\n',sep='')
                cat('     Scaled max. abs. change in coef.: ',
                    format(round(scrit[3],digits)),'\n',sep='')
                cat('     Sum abs. val. of grad. vector.: ',
                    format(round(scrit[4],digits)),'\n',sep='')
        }

        if(scrit[icrit] < tolerance) {
                        converged <- TRUE
                        nstep <- sd.step
                        break
        }
        if(sd.step >= itmax) {
                cat(paste0("Failed to converge in ",itmax,
                          " steepest descent steps.\n",sep=""))
                converged <- FALSE
                nstep     <- sd.step
                break
        }
        oath    <- theta
        old.ll  <- ll
        sd.step <- sd.step + 1
}

tpm   <- getTpm(theta,K,stationary=TRUE)
ispd  <- revise.ispd(tpm)
Rho   <- getRho(theta,K,rhovals=rhovals,stationary=TRUE,prednames="Intercept")
xxx   <- get.hgl(theta,K,y,hglmethod=hglmethod)
hess  <- xxx$hess
grad  <- xxx$grad
ll    <- xxx$ll
AIC   <- -2*ll+2*npar
BIC   <- -2*ll+bicm*npar
rownames(hess) <- colnames(hess) <- names(theta)

names(scrit) <- c("PCLL","L2","Linf","ABSGRD")
list(Rho=Rho,tpm=tpm,ispd=ispd,log.like=ll,hessian=hess,grad=grad,
     stopCrit=scrit,par0=par0,npar=npar,bicm=bicm,converged=converged,
     nstep=nstep,AIC=AIC,BIC=BIC)
}
