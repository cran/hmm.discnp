revise.tpm <- function(xi,mixture) {
	if(mixture)  matrix(apply(xi,2,sum)/sum(xi),byrow=T,nrow=nrow(xi))
	else xi/apply(xi,1,sum)
}