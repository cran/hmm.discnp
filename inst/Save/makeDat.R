makeDat <- function(y,X,newstyle,lvls) {
    if(newstyle) {
        if(is.null(X)) {
            aI <- TRUE
            pnms <- "Intercept"
        } else {
            pnms <- attr(X,"prednames")
            aI   <- pnms[[1]] == "Intercept"
        }
        n    <- length(y)
        rslt <- vector("list",n)
        for(i in seq(along=rslt)) {
            yi <- factor(y[[i]],levels=lvls)
            Xi <- if(aI) cbind(1,X[[i]]) else X[[i]]
            rslt[[i]] <- data.frame(yi,Xi)
            names(rslt[[i]]) <- c("y",pnms)
        }
        return(rslt)
    }
    return(lapply(y,function(x,lvls){factor(x,levels=lvls)},lvls=lvls))
} 
