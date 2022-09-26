steepest <- function(K,y,theta) {
foo <- function(x,K,y,theta,gv) { 
		get.l(theta+x*gv,K,y)
        } 
gv     <- get.gl(theta,K,y)$grad
gv     <- gv/sqrt(sum(gv^2))
fooMax <- optimize(foo,c(0,1),maximum=TRUE,K=K,y=y,theta=theta,gv=gv)
# Numerical inaccuracy can result in a "maximum" that is
# a miniscule amount *less* than than the log likelihood value
# at the original theta.  (In such circumstances the original
# theta is a local maximum.)  Using the maximum thus produced
# by optimize() would result in a (small) *decrease* in the
# log likelihood.  We guard against this as follows:
con   <- fooMax$maximum
lwbnd <- foo(0,K,y,theta,gv)
if(fooMax$objective < lwbnd) con <- 0
theta + con*gv
}
