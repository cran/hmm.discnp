sim.hmm <- function(nsim,tpm,Rho,nrep=1) {
#
# Function sim.hmm to simulate data from a hidden Markov
# model with transition probability matrix tpm, and discrete
# (non-parametric) distributions specified by the matrix Rho.
#

# Check for validity of tpm and Rho arguments:
if(ncol(tpm) != nrow(tpm))
	stop("The matrix tpm must be square.\n")
if(ncol(tpm) != ncol(Rho))
	stop("Mismatch in dimensions of tpm and Rho.\n")
if(any(tpm<0)) stop("Negative entries in tpm.\n")
if(any(Rho<0)) stop("Negative entries in Rho.\n")
xxx <- apply(tpm,1,sum)
if(!identical(all.equal(xxx,rep(1,nrow(tpm))),TRUE))
	stop("Rows of tpm do not all sum to 1.\n")
xxx <- apply(Rho,2,sum)
if(!identical(all.equal(xxx,rep(1,ncol(Rho))),TRUE))
	stop("Columns of Rho do not all sum to 1.\n")

ispd <- revise.ispd(tpm)
K    <- ncol(Rho)
M    <- nrow(Rho)
rslt <- list()
for(j in 1:nrep) {
	s1     <- sample(1:K,1,prob=ispd)
	y      <- list()
	y[[1]] <- sample(1:M,1,prob=Rho[,s1])
	for(i in 2:nsim) {
		s1 <- sample(1:K,1,prob=tpm[s1,])
		y[[i]] <- sample(1:M,1,prob=Rho[,s1])
	}
	rslt[[j]] <- unlist(y)
}
if(nrep>1) matrix(unlist(rslt),ncol=nrep) else unlist(rslt)
}
