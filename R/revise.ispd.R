revise.ispd <- function(tpm)
{
# Function revise.ispd.  To revise the initial state probability
# distribution in the setting wherein ispd is taken to be the steady
# state distribution.

eee   <- eigen(t(tpm))
k <- match(1,round(eee$values,6))
if(length(k) != 1) {
	cat('Problems with eigenvalues:\n')
	print(eee$values)
	stop()
}
ispd <- Re(eee$vectors[,k])
ispd/sum(ispd)
}
