fix.tpm <- function(x,K) {
y <- matrix(x,nrow=K,byrow=TRUE)
z <- exp(cbind(y,0))
z/apply(z,1,sum)
}
