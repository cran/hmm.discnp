sim.hmm <- function(nsim,tpm,Rho,nrep=1) {
#
# Function sim.hmm to simulate data from a hidden Markov
# model with transition probability matrix tpm, and discrete
# (non-parametric) distributions specified by the matrix Rho.
#

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
