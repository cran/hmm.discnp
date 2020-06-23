viterbi <- function(y,model=NULL,tpm=NULL,Rho=NULL,ispd=NULL,
                    log=FALSE,warn=TRUE) {
#
# Function viterbi to apply the Viterbi algorithm to a collection
# of data sequences, given the parameters of the model.
#

if(inherits(y,"multipleHmmDataSets")) {
    rslt <- lapply(y,function(x,model,tpm,Rho,ispd,log,warn){
                         viterbi(x,model,tpm,Rho,ispd,log,warn)},
                         model=model,tpm=tpm,Rho=Rho,ispd=ispd,
                         log=log,warn=warn)
    return(rslt)
}

# If ``model'' is present, get the parameters from that, and
# ignore those (if any) specified as separate arguments.
if(!is.null(model)) {
	tpm  <- model$tpm
	Rho  <- model$Rho
	ispd <- model$ispd
}
if(is.null(tpm) | is.null(Rho))
    stop("At least one of \"tpm\" and \"Rho\" was not supplied.\n")

# Convert Rho if necessary.
if(inherits(Rho,"matrix")) Rho <- cnvrtRho(Rho)

# Set the type:
# 1 <--> univariate
# 2 <--> bivariate, independent
# 3 <--> bivariate, dependent
if(inherits(Rho,"data.frame")) {
    type <- 1
} else if(inherits(Rho,"list")) {
    type <- 2
} else if(inherits(Rho,"array")) {
    type <- 3
} else {
    stop("Object \"Rho\" has an incorrect class.\n")
}
K <- nrow(tpm)
if(missing(y)) {
	y <- if(!is.null(model)) model[["y"]] else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
}

# Make sure y is a list, and get the number of sequences and
# lengths of these sequences.
y    <- tidyList(y)
nseq <- length(y)
lns  <- sapply(y,length)
lvls <- attr(y,"lvls")

# Build ispd if it was given as NULL
if(is.null(ispd)) ispd <- revise.ispd(tpm)

# Make sure that the y-values are compatible with Rho.
Rho <- check.yval(attr(y,"lvls"),Rho,type,warn=warn)

fys <- function(y,Rho,type,lvls) {
   Dat <- switch(EXPR=type,
               list(data.frame(y=factor(y,levels=lvls),Intercept=1)),
               list(data.frame(y1=factor(y[1],levels=lvls[[1]]),
                               y2=factor(y[2],levels=lvls[[2]]),
                               Intercept=1)),
               list(data.frame(y1=factor(y[1],levels=lvls[[1]]),
                               y2=factor(y[2],levels=lvls[[2]]),
                               Intercept=1)))
   as.vector(ffun(Dat,Rho,type))
}

# Set the names of the states.
stnms <- switch(type,
             levels(Rho$state),
             colnames(Rho[[1]]),
             dimnames(Rho)[[3]]
         )
if(is.null(stnms)) stnms <- as.character(1:nrow(tpm))

# Rubber hits road.
rslt <- vector("list",nseq)
for(j in 1:nseq) {
	psi <- vector("list",nseq)
        yt  <- y[[j]][1,]
        rrr <- fys(yt,Rho,type,lvls)
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
                rrr <- fys(yt,Rho,type,lvls)
		if(log) {
                    delta <- log(rrr) + apply(delta + log(tpm),2,max)
                } else {
		    delta <- rrr*apply(delta*tpm,2,max)
                    delta <- delta/sum(delta)
                }
	}
	temp <- vector("list",nj)
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
        sss <- stnms[unlist(temp[[1]])]
        sss <- matrix(sss, nrow = nj)
        rslt[[j]] <- if(ncol(sss)==1) as.vector(sss) else sss
}
if(nseq==1) rslt[[1]] else rslt
}
