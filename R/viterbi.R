viterbi <- function(y,model=NULL,tpm,Rho,ispd=NULL,log=FALSE,warn=TRUE) {
#
# Function viterbi to apply the Viterbi algorithm to a collection
# of data sequences, given the parameters of the model.
#

# If ``model'' is present, get the parameters from that, and
# ignore those (if any) specified as separate arguments.
if(!is.null(model)) {
	tpm  <- model$tpm
	Rho  <- model$Rho
	ispd <- model$ispd
}

# Set the type:
type <- switch(class(Rho),matrix=1,list=2,array=3,NULL)
if(is.null(type)) stop("Argument \"Rho\" is not of an appropriate form.\n")

K <- nrow(tpm)
if(missing(y)) {
	y <- if(!is.null(model)) model$y else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
}

# Make sure y is a list, and get the number of sequences and
# lengths of these sequences.
y    <- tidyList(y)
nseq <- length(y)
lns  <- sapply(y,length)

# Build ispd if it was given as NULL
if(is.null(ispd)) ispd <- revise.ispd(tpm)

# Make sure that the y-values are compatible with Rho.
Rho <- check.yval(y,Rho,type,warn=warn)

fys <- function(y,s,Rho,type) {
    switch(type,Rho[y,s],Rho[[1]][y[1],s]*Rho[[2]][y[2],s],
                Rho[cbind(y[1],y[2],s)])
}

sK <- switch(type, colnames(Rho),colnames(Rho[[1]]),dimnames(Rho)[[3]])
if(is.null(sK)) sK <- 1:K
rslt <- list()
for(j in 1:nseq) {
	psi <- vector("list",nseq)
        yt  <- y[[j]][1,]
        rrr <- fys(yt,sK,Rho,type)
        if(log) {
            delta <- log(ispd) + log(rrr)
        } else {
	    delta <- ispd*rrr
	    delta <- delta/sum(delta)
        }
	nj <- lns[j]
	for(tt in 2:nj) {
		if(log) {
                    tmp <- apply(delta + log(tpm),2,
                             function(x){((1:length(x))[x==max(x)])}
                             )
                } else {
		    tmp <- apply(delta*tpm,2,
                             function(x){((1:length(x))[x==max(x)])}
                             )
                }
	        psi[[tt]] <- tmp # Note that tmp will be a list of
		                 # vectors, each of length between
                                 # 1 and K = the number of states.
                yt  <- y[[j]][tt,]
                rrr <- fys(yt,sK,Rho,type)
		if(log) {
                    delta <- log(rrr) + apply(delta + log(tpm),2,max)
                } else {
		    delta <- rrr*apply(delta*tpm,2,max)
                    delta <- delta/sum(delta)
                }
	}
	temp <- list()
	temp[[nj]] <- (1:K)[delta==max(delta)]
	for(tt in (nj-1):1) {
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
        rrr <- matrix(unlist(temp[[1]]), nrow = nj)
        rslt[[j]] <- if(ncol(rrr)==1) as.vector(rrr) else rrr
}
if(nseq==1) rslt[[1]] else rslt
}
