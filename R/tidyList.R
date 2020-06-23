tidyList <- local({

goodsort <- function(x){
   x <- unique(x)
   x <- x[!is.na(x)]
   if(is.numeric(x)) return(sort(x)) # Can't happen now 30/06/2018.
   xn <- suppressWarnings(as.numeric(x))
   if(!any(is.na(xn))) return(x[order(xn)])
   sort(x)
}

function(y,rp=c("response","predictor"),addIntercept=NULL,yval=NULL) {

# If "y" is already a "tidyList" do nothing.
if(inherits(y,"tidyList")) return(y)

if(!inherits(y,"list")) y <- list(y)
rp <- match.arg(rp)

ArgNm <- if(rp=="response") "y" else "X"
ddd   <- lapply(y,dim)
nnn   <- sapply(ddd,is.null)
if(any(nnn)) {
    if(all(nnn)) {
        y <- lapply(y,matrix)
# Note that if any of the entries of y were factors, the
# foregoing line will turn them into (one-column) character
# matrices.
    } else {
        stop(paste0("Argument \"",ArgNm,"\" is a list of objects",
                    " with differing structure.\n"))
    }
}

# The object y is now a list of objects all of which have dimensions.
if(!all(sapply(y,function(x){inherits(x,c("matrix","data.frame"))})))
    stop(paste0("At least one entry of \"",ArgNm,"\" is neither",
                " a matrix nor a data frame.\n"))
ncy <- unique(sapply(y,ncol))
if(length(ncy) > 1) stop(paste0("Entries of \"",ArgNm,"\" have differing",
                                " numbers of columns.\n"))

y <- lapply(y,as.matrix)
# Note that the entries in the y list were matrices or
# or data frames.  If any of the columns of the data frames
# were factors, the foreging line will coerce these data frames
# into *character* matrices.

if(rp=="predictor") {
# Get the predictor names from the column names of the
# entries of y.
    colnms    <- lapply(y,colnames)
    nullNames <- sapply(colnms,is.null)
    if(all(nullNames)) {
        prednames <- NULL
    } else {
        if(!any(nullNames)) {
            mcnms <- matrix(unlist(colnms),byrow=TRUE,ncol=ncy)
            mcnms <- mcnms[!duplicated(mcnms),,drop=FALSE]
            if(nrow(mcnms)==1) {
                prednames <- mcnms[1,]
                ok        <- TRUE
            } else ok <- FALSE
        } else ok <- FALSE
        if(!ok) stop("All predictor matrices must have the same column names.\n")
    }

    hasInt <- checkConstColumns(y,prednames)

    if(is.null(prednames)) {
        lout <- if(hasInt) ncy-1 else ncy
        iseq <- seq(from=1,length.out=lout)
        prednames <- if(length(iseq)) paste0("V",iseq) else NULL
        if(hasInt) prednames <- c("Intercept",prednames)
    }

    if(is.null(addIntercept)) addIntercept <- !hasInt
    if(addIntercept) {
        if(hasInt) {
            whinge <- paste0("Argument \"addIntercept\" is TRUE, but the\n",
                             " predictors already have an intercept.\n")
            stop(whinge)
        } else {
            prednames <- c("Intercept",prednames)
        }
        y <- lapply(y,function(x){cbind(1,x)})
    }
    y <- lapply(y,function(x,prednames){colnames(x) <- prednames;x},
                  prednames=prednames)
    ncy <- ncy+addIntercept
    attr(y,"prednames") <- prednames
    attr(y,"ncol") <- ncy
    class(y) <- c("tidyList","list")
    return(y)
}

# We're now just looking at the response.  Set the "ncol" attribute.
# (This will be either 1 or 2.)
attr(y,"ncol") <- ncy

# Note that we need to coerce the y-values (the responses) to
# character mode, otherwise we may/will wind up with mismatches to
# the dimension names of Rho.
my <- unique(sapply(y,mode))

# At this stage any entry of y is a matrix of mode either
# "numeric" or "character".  If they are all numeric, then
# we can say that the original data were numeric.  If they
# are all character, then we can't.  If there is a mix of the
# two, then alles upgefucken ist.
if(length(my) > 1) stop("Mixture of modes in data.\n")
odn <- my=="numeric" # "Original data numeric".
y   <- lapply(y,function(x){mode(x) <- "character"; x})

# Univariate.
if(ncy==1) {
    obsVal <- goodsort(unlist(y))
    if(is.null(yval)) {
        lvls <- obsVal
    } else {
        if(!all(obsVal%in%yval))
            stop("Specified y values do not include all observed y values.\n")
        lvls <- as.character(yval)
    }
    attr(y,"lvls") <- lvls
    attr(y,"parity") <- "univar"
    attr(y,"numeric") <- odn
    class(y) <- c("tidyList","list")
    return(y)
}

# Bivariate.
if(ncy==2) {
    obsVal1 <- goodsort(unlist(lapply(y,function(x){x[,1]})))
    obsVal2 <- goodsort(unlist(lapply(y,function(x){x[,2]})))
    if(is.null(yval)) {
        lvls  <- list(obsVal1,obsVal2)
    } else {
        if(!(is.list(yval) && length(yval)==2))
                stop("Argument \"yval\" is not of the right form.\n")
        if(!(all(obsVal1%in%yval[[1]]) & all(obsVal2%in%yval[[2]])))
            stop("Specified y values do not include all observed y values.\n")
        lvls <- lapply(yval,as.character)
    }
    attr(y,"lvls") <- lvls
    attr(y,"parity") <- "bivar"
    attr(y,"numeric") <- odn
    class(y) <- c("tidyList","list")
    return(y)
}

whinge <- paste("Argument \"y\" consists of matrices with more",
                "than two columns.\n")
stop(whinge)
}})
