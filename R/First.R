.First.lib <- function(lib,pkg) {
	library.dynam("hmm.discnp", pkg, lib)
	ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
	cat(paste(pkg, ver, "\n\n"))
	cat("WARNING:  The package has changed substantially\n")
	cat("from the previous version.  New functions have been\n")
	cat("added and both the argument lists and the returned\n")
	cat("values from old functions have new forms.  Please\n")
	cat("read the ChangeLog and the documentation.\n\n")
}
