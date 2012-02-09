check.yval <- function(y,Rho) {
yval  <- unique(unlist(y))
fname <- as.character(sys.call(-1))[1]
if(is.na(fname)) fname <- "call from the command line"
if(is.null(row.names(Rho))) {
	yval <- as.numeric(yval)
	OK   <- all(yval%in%(1:nrow(Rho)))
	if(!OK) stop(paste("In ",fname,".  The values of \"y\" must be in ",
                     "\"1:nrow(Rho)\".\n",sep=""),call.=FALSE)
} else {
	yval <- as.character(yval)
	OK   <- all(yval%in%row.names(Rho))
	if(!OK) stop(paste("In ",fname,".  The values of \"y\" must be in ",
                     "\"row.names(Rho)\".\n",sep=""),call.=FALSE)
}
return(invisible())
}
