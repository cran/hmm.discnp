.onLoad <- function(lib, pkg) {
	library.dynam("hmm.discnp", pkg, lib)
	ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
        packageStartupMessage(paste(pkg, ver))
	packageStartupMessage(
            paste("\n     PLEASE NOTE:  The package has changed substantially",
                  "\n     from the 0.0-x versions.  New functions have been",
                  "\n     added and both the argument lists and the returned",
                  "\n     values from old functions have new forms.  Please",
                  "\n     read the ChangeLog and the documentation.\n")
        )
}
