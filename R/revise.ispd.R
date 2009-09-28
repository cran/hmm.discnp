revise.ispd <- function(tpm=NULL,gamma=NULL,lns=NULL) {
# Function revise.ispd.  To revise the initial state probability
# distribution

if(is.null(tpm) + (is.null(gamma) & is.null(lns)) != 1)
	stop(paste("Either \"tpm\" should be given OR\n",
                   "\"gamma\" and \"lns\" should be given.\n"))

# If ispd is taken to be the steady state distribution:
if(!is.null(tpm)) {
	eee   <- eigen(t(tpm))
	k <- match(1,round(eee$values,6))
	if(length(k) != 1) {
		cat('Problems with eigenvalues:\n')
		print(eee$values)
		stop()
	}
	v <- Re(eee$vectors[,k])
	return(v/sum(v))
}

# Steady state not assumed:
v <- 0
j <- 1
nseq <- length(lns)
for(i in 1:nseq) {
		v <- v + gamma[,j]
		j <- j + lns[i]
}
v/nseq 
}
