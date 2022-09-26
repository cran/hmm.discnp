lmstep <- function(theta,K,y,lmc,hglmethod) {
#
# Function lmstep to do one Levenberg-Marquardt step.
#
xxx     <- get.hgl(theta,K,y,hglmethod=hglmethod)
hess    <- xxx$hess
ev      <- eigen(hess)$values
osag    <- sum(abs(xxx$grad))
npar    <- length(theta)
nms     <- names(theta)
old.ll  <- xxx$ll
steepit <- FALSE
repeat {
	lll       <- lmc + max(0,ev)
	hess.mod  <- hess - lll*diag(npar)
        invhm     <- try(solve(hess.mod),TRUE)
        if(inherits(invhm,"try-error")) {
            steepit <- TRUE
            break
        }
	theta.new <- theta - invhm%*%xxx$grad
	yyy       <- try(get.gl(theta.new,K,y),TRUE)
	if(inherits(yyy,"try-error")) {
		steepit <- TRUE
		break
	}
	new.ll    <- yyy$ll
        nsag      <- sum(abs(yyy$grad))
        if(all(is.finite(c(old.ll,new.ll,osag,nsag))) &&
                (new.ll > old.ll & nsag < osag)) break
	lmc <- 10*lmc
        if(1/lmc < sqrt(.Machine$double.eps)) {
		steepit <- TRUE
		break
	}
}
if(steepit) {
	theta.new <- steepest(K,y,theta)
	yyy <- get.gl(theta.new, K, y)
	lmc <- 1
}
theta.new <- as.vector(theta.new)
names(theta.new) <- nms
grad <- yyy$grad
names(grad) <- nms
list(theta=theta.new,ll=yyy$ll,grad=grad,lmc=lmc,used.steepest=steepit)
}
