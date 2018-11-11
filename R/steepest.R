steepest <- function(K,y,theta){
foo <- function(x,K,y,theta,gv) { 
		get.l(theta+x*gv, K, y)
        } 
gv <- get.gl(theta,K,y)$grad
gv <- gv/sqrt(sum(gv^2))
poop <- optimize(foo,c(-1,1),maximum=TRUE,K=K,y=y,
                 theta=theta,gv=gv)
con <- poop$maximum
theta + con*gv
}
