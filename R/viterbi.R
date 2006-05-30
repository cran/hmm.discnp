viterbi <- function(y,object=NULL,tpm,Rho,ispd,yval=NULL) {
#
# Function vitervi to apply the Viterbi algorithm to a collection
# of equi-length data sequences, given the parameters of the model.
#

y  <- tidyup(y,yval)
n  <- nrow(y)
nc <- ncol(y)

if(!is.null(object)) {
	tpm  <- object$tpm
	Rho  <- object$Rho
	ispd <- object$ispd
}
K <- nrow(tpm)

rslt <- list()
for(j in 1:nc) {
	psi <- list()
	delta <- ispd*Rho[y[1,j],]
	delta <- delta/sum(delta)
	for(tt in 2:n) {
		tmp <- apply(delta*tpm,2,
                             function(x){((1:length(x))[x==max(x)])}
                             )
		psi[[tt]] <- tmp # Note that tmp will be a list of
		                 # vectors, each of length between
                                 # 1 and K = the number of states.
		delta <- Rho[y[tt,j],]*apply(delta*tpm,2,max)
	}

	temp <- list()
	temp[[n]] <- (1:K)[delta==max(delta)]
	for(tt in (n-1):1) {
		i <- 0
		temp[[tt]] <- list()
		for(x in temp[[tt+1]]) {
			k <- x[1]
			for(w in psi[[tt+1]][[k]]) {
				i <- i+1
				temp[[tt]][[i]] <- c(w,x)
			}
		}
	}
        rrr <- matrix(unlist(temp[[1]]), nrow = n)
        rslt[[j]] <- if(ncol(rrr)==1) as.vector(rrr) else rrr
}
if(nc==1) rslt <- rslt[[1]]
rslt
}
