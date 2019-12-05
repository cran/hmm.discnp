mps <- function(y,model=NULL,tpm,Rho,ispd=NULL,warn=TRUE){
#
# Function mps: most probable states.
#
# To do:  Something about supplying newstyle and X. (06/04/2017)
#

if(!is.null(model)) {
	tpm  <- model$tpm
	Rho  <- model$Rho
	ispd <- model$ispd
}

# If Rho is presented in the "logistic style" parameterisation,
# convert it to a matrix of probabilities.
if(inherits(Rho,"data.frame")) Rho <- cnvrtRho(Rho)

stnms <- rownames(tpm)
if(is.null(ispd)) ispd <- revise.ispd(tpm)
if(missing(y)) {
	y <- if(!is.null(model)) model$y else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
}
y <- tidyList(y)
y <- makeDat(y,X=NULL)

# Set the type:
if(inherits(Rho,"data.frame")) {
    type <- 1
} else if(inherits(Rho,"list")) {
    type <- 3
} else if(inherits(Rho,c("matrix","array"))) {
    if(length(dim(Rho))==2) type <- 2
    else if(length(dim(Rho))==3) type <- 4
    else stop("Object \"Rho\" can be of dimension 2 or 3 only.\n")
} else {
    stop("Object \"Rho\" has an incorrect class.\n")
}
# Note: type = 1 can't happen here; Rho would have been converted
# from data frame to matrix.

Rho  <- check.yval(attr(y,"lvls"),Rho,type,warn=warn)
lns  <- sapply(y,nrow)
nseq <- length(y)
fy   <- ffun(y,Rho,type)
rp   <- recurse(fy, tpm, ispd, lns)
xxx  <- apply(rp$gamma, 2, which.max)
if(!is.null(stnms)) xxx <- stnms[xxx]
if(nseq==1) return(xxx)
rslt  <- list()
jstop <- 0
for(i in 1:nseq) {
        jstart <- jstop+1
        jstop  <- jstop + lns[i]
	rslt[[i]] <- xxx[jstart:jstop]
}
rslt
}
