.First.lib <- function(lib,pkg) {
	library.dynam("hmm.discnp", pkg, lib)
	ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
	cat(paste(pkg, ver, "\n"))
}
