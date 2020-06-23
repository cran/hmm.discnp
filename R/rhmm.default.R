rhmm.default <- local({

canBeNumeric <- function(x) {
    !any(is.na(suppressWarnings(as.numeric(x))))
}

function(model,...,nsim=1,verbose=FALSE,ylengths,
                         nafrac=NULL,fep=NULL,tpm,Rho,
                         ispd=NULL,yval=NULL,drop=TRUE,forceNumeric=TRUE) {
#
# Function rhmm.default to simulate data from a hidden Markov
# model with transition probability matrix tpm, and discrete
# (non-parametric) distributions specified by Rho which
# may be a matrix (univariate emissions), a list of two
# matrices (independent bivariate emissions) or an array
# (dependent bivariate emissions).
#

# Check on compatibility of ylengths and ispd:
if(missing(ylengths)) stop("Argument \"ylengths\" must be supplied.\n")
nseq <- length(ylengths)
if(inherits(ispd,"matrix")) {
   if(nseq != ncol(ispd)) {
       whinge <- paste("Number of columns of \"ispd\" is not equal to",
                       "the length of \"ylengths\".\n")
       stop(whinge)
   }
}

# Convert Rho if necessary.
if(inherits(Rho,"matrix")) Rho <- cnvrtRho(Rho)

# Set the type:
# 1 = univariate, 2 = bivariate independent, 3 = bivariate dependent.
if(inherits(Rho,"data.frame")) {
   type <- 1
} else if(inherits(Rho,"list")) {
   type <- 2
} else if(inherits(Rho,"array")) {
   type <- 3
} else {
    stop("Argument \"Rho\" is not of an appropriate form.\n")
}
parity <- if(type==1) "univar" else "bivar"
nyv    <- is.null(yval)
if(nyv & type > 1) yval <- vector("list",2)

# Check for validity of the Rho argument and reconcile dimension
# names with "yval".
switch(type,
    {
        if(is.null(yval)) {
            yval <- levels(Rho$y)
        } else {
            if(length(yval)!=length(levels(Rho$y)))
                stop(paste("Mismatch between length of \"yval\"",
                           "and number of y values specified by \"Rho\".\n"))
        }
        Rho <- cnvrtRho(Rho)
        rownames(Rho) <- yval
    },
    {
      if(ncol(Rho[[1]]) != ncol(Rho[[2]]))
          stop("Mismatch in number of states between Rho[[1]] and Rho[[2]].\n")
      if(length(yval) !=2) # If yval was NULL it would have been set to
                           # an empty list of length 2.
          stop("Argument \"yval\", if not NULL, should be a list of length 2.\n")
      for(j in 1:2) {
        if(any(Rho[[j]] < 0))
            stop(paste("Negative entries in Rho[[",j,"]].\n"),sep="")
        xxx <- unname(apply(Rho[[j]],2,sum))
        if(!isTRUE(all.equal(xxx,rep(1,ncol(Rho[[j]])))))
            stop(paste("Columns of Rho[[",j,"]] do not all sum to 1.\n"),sep="")
        if(nyv) {
           yval[[j]] <- if(is.null(rownames(Rho[[j]]))) 1:nrow(Rho[[j]])
                             else rownames(Rho[[j]])
        } else {
	  if(length(yval[[j]])!=nrow(Rho[[j]]))
		  stop(paste("Mismatch between length of yval[[",j,"]]",
                             "and number of rows of Rho[[",j,"]].\n"),sep="")
       }
       rownames(Rho[[j]]) <- yval[[j]]
     }
    },
    {
      if(any(Rho<0)) stop("Negative entries in Rho.\n")
      xxx <- unname(apply(Rho,3,sum))
      if(!isTRUE(all.equal(xxx,rep(1,dim(Rho)[3]))))
	stop("Layers of Rho do not all sum to 1.\n")
      for(j in 1:2) {
        if(nyv) {
           yval[[j]] <- if(is.null(dimnames(Rho)[[j]])) 1:nrow(Rho[[j]])
                             else dimnames(Rho)[[j]]
        } else {
	  if(length(yval[[j]])!=dim(Rho)[j])
              stop(paste("Mismatch between the length of yval[[",j,"]]",
                         "and the ",j,"th dimension of Rho.\n"),sep="")
       }
      }
    })


K <- switch(type,ncol(Rho),ncol(Rho[[1]]),dim(Rho)[3])

if(inherits(yval,"list")) {
    ok <- sapply(yval,canBeNumeric)
} else ok <- canBeNumeric(yval)

if(length(ok) == 1) {
    obsRnum <- ok
} else if(length(ok) == 2) {
    nok <- sum(ok)
    if(nok == 2) {
        obsRnum <- TRUE
    } else if(nok == 0) {
        obsRnum <- FALSE
    } else stop("Inconsistent modes for y-values.\n")
}
yval.save <- yval
if(obsRnum & forceNumeric) {
    if(inherits(yval,"list")) {
        yval <- lapply(yval,as.numeric)
    } else {
        yval <- as.numeric(yval)
    }
}

xample <- function(n,yval,Rho,state,type) {
# We could calculate "type" here but it's faster (???)
# to use the pre-calculated value.
switch(type,
    return(matrix(sample(yval,size=n,prob=Rho[,state],replace=TRUE))),
    {
        y1 <- sample(yval[[1]],size=n,prob=Rho[[1]][,state],replace=TRUE)
        y2 <- sample(yval[[2]],size=n,prob=Rho[[2]][,state],replace=TRUE)
        return(cbind(y1,y2))
    },
    {
    mn <- prod(dim(Rho)[1:2])
    ij <- sample(1:mn,size=n,prob=as.vector(Rho[,,state]), replace=TRUE)
    i  <- row(Rho[,,state])[ij]
    j  <- col(Rho[,,state])[ij]
    return(cbind(yval[[1]][i],yval[[2]][j]))
    })
}

rslt <- vector("list",nsim)
# If there is only one state generate i.i.d. data.
if(K == 1) {
    for(i in 1:nsim) {
        rslt[[i]] <- lapply(ylengths,xample,yval=yval,Rho=Rho,state=1,type=type)
    }
    return(if(drop & nsim==1) rslt[[1]] else rslt)
}

# More than one state; check on tpm:
if(ncol(tpm) != nrow(tpm))
	stop("The matrix tpm must be square.\n")
if(ncol(tpm) != K)
	stop("Mismatch between dimensions of tpm and Rho.\n")
if(any(tpm<0)) stop("Negative entries in tpm.\n")
xxx <- unname(apply(tpm,1,sum))
if(!isTRUE(all.equal(xxx,rep(1,nrow(tpm)))))
	stop("Rows of tpm do not all sum to 1.\n")

if(is.null(ispd)) ispd <- revise.ispd(tpm)
jr   <- 0
for(i in 1:nsim) {
    tempRes <- vector("list",nseq)
    for(j in 1:nseq) {
        jr    <- jr+1
        M     <- matrix(if(obsRnum) 0 else "",nrow=ylengths[j],
                        ncol=if(type==1) 1 else 2)
        prb   <- if(inherits(ispd,"matrix")) ispd[,j] else ispd
	s1    <- sample(1:K,1,prob=prb)
        M[1,] <- xample(1,yval,Rho,state=s1,type)
	for(k in 2:ylengths[j]) {
            jr <- jr+1
            s1 <- sample(1:K,1,prob=tpm[s1,])
            M[k,] <- xample(1,yval,Rho,state=s1,type)
            if(verbose) {
                if(jr%%1000 == 0) cat(jr,"")
                if(jr%%10000 == 0) cat("\n")
            }
	}
        tempRes[[j]] <- M
    }
    if(!is.null(nafrac)) tempRes <- misstify(tempRes,nafrac=nafrac,fep=fep)
    if(type==1) tempRes <- lapply(tempRes,as.vector)
    if(nseq==1 & drop) tempRes <- tempRes[[1]]
    attr(tempRes,"uval") <- yval.save
    attr(tempRes,"parity") <- parity
    rslt[[i]] <- tempRes
}
if(verbose) cat("\n")
if(nsim==1 & drop) return(rslt[[1]])
class(rslt) <- c(class(rslt),"multipleHmmDataSets")
rslt
}
})
