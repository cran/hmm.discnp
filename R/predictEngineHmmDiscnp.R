predictEngineHmmDiscnp <- function(stateProbs,Rho,numb,drop) {

# Convert Rho if necessary (just for consistency with the
# overall pattern; we'll actually convert it *back* again later).
if(inherits(Rho,"matrix")) Rho <- cnvrtRho(Rho)

# Set the type:
if(inherits(Rho,"data.frame")) {
    type <- 1
} else if(inherits(Rho,"list")) {
    type <- 2
} else if(inherits(Rho,"array")) {
    type <- 3
} else {
    stop("Object \"Rho\" has an incorrect class.\n")
}

if(!inherits(stateProbs,"list")) stateProbs <- list(stateProbs)
fitVal <- switch(type,
# Univariate
    {
        Roe <- cnvrtRho(Rho)
        if(numb) {
            yval <- as.numeric(levels(Rho$y))
            if (any(is.na(yval))) {
                whinge <- paste0("Argument \"object\" indicates that",
                                 " the y-values are numeric,\n",
                                 "but the levels of Rho$y are not numeric.\n")
                stop(whinge)
            } 
            cmns <- apply(yval * Roe, 2, sum)
            lapply(stateProbs,function(x,cmns){apply(cmns*x,2,sum)},
                             cmns=cmns)
        } else {
            lapply(stateProbs,function(x,Rho){t(Rho%*%x)},Rho=Roe)
        }
    },

# Bivariate independent:
    {
        if(numb) {
            lapply(stateProbs,function(x,Rho) {
                temp <- vector("list",2)
                for (j in 1:2) {
                    yval <- as.numeric(row.names(Rho[[j]]))
                    if (any(is.na(yval))) {
                        whinge <- paste0("Argument \"object\" indicates that",
                                         " the y-values are numeric,\n",
                                         "but the row names of Rho[[",j,"]]",
                                         " are not numeric.\n")
                        stop(whinge)
                    } 
                    cmns <- apply(yval * Rho[[j]], 2, sum)
                    temp[[j]] <- apply(cmns*x,2,sum)
                }
                cbind(temp[[1]],temp[[2]])},Rho=Rho)
        } else {
            lapply(stateProbs,function(x,Rho){
                temp <- vector("list",2)
                for(j in 1:2) {
                    temp[[j]] <- t(Rho[[j]]%*%x)
                }
                aaa <- array(NA,c(nrow(Rho[[1]]),nrow(Rho[[2]]),ncol(x)))
                aaa[1,,] <- temp[[1]]
                aaa[,2,] <- temp[[2]]
                aaa},Rho=Rho)
        }
    },

# Bivariate dependent:
    {
        if(numb) {
            lapply(stateProbs,function(x,Rho){
                yval <- vector("list", 2)
                for(j in 1:2) {
                    yval[[j]] <- as.numeric(dimnames(Rho)[[j]])
                    if(any(is.na(yval[[j]]))) { 
                        whinge <- paste0("Argument \"object\" indicates that",
                                         " the y-values are numeric,\n",
                                         "but dimnames(Rho)[",j,"]",
                                         " is not numeric.\n")
                        stop(whinge)
                    }
                } 
                temp <- vector("list",2)
                for (j in 1:2) {
                    RT <- apply(Rho, c(j, 3), sum)
                    cmns <- apply(yval[[j]] * RT, 2, sum)
                    temp[[j]] <- apply(cmns * x, 2, sum)
                }
                cbind(temp[[1]],temp[[2]])},Rho=Rho)
        } else {
            lapply(stateProbs,function(x,Rho){
                K   <- dim(Rho)[3]
                aaa <- vector("list",K)
                for(k in 1:K) {
                    aaa[[k]] <- outer(Rho[,,k],x[k,])
                }
                Reduce("+",aaa)},Rho=Rho) 
        }
   })
   if(length(fitVal)==1 & drop) fitVal <- fitVal[[1]]
   fitVal
}
