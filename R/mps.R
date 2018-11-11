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
stnms <- rownames(tpm)
if(is.null(ispd)) ispd <- revise.ispd(tpm)
if(missing(y)) {
	y <- if(!is.null(model)) model$y else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
}
y    <- tidyList(y)
type <- if(is.matrix(Rho)) {
            1
        } else if(is.list(Rho)) {
            2
        } else if(is.array(Rho)) 3

if(is.null(type)) stop("Argument \"Rho\" is not of an appropriate form.\n")
Rho  <- check.yval(y,Rho,type,warn=warn)
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
