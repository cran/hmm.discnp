revise.rho <- function(y,gamma,m) {
	there <- !is.na(c(y))
	t1 <- apply(gamma[,there],1,
		function(x,index){tapply(x,index,sum)},c(y)[there])
	t2 <- t(t(t1)/apply(t1,2,sum))
	rslt <- matrix(0,nrow=m,ncol=nrow(gamma))
	rslt[sort(unique(as.vector(y))),] <- t2
	rslt
}
