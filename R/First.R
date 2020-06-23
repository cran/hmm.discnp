.onLoad <- function(lib, pkg) {
	library.dynam("hmm.discnp", pkg, lib)
}

.onAttach <- function(lib, pkg) {
	ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
        packageStartupMessage(paste(pkg, ver))
	msg <- paste("\n     This package has changed SUBSTANTIALLY from its",
                     "\n     previous release.  Read the documentation",
                     "\n     carefully.  Note in particular that there is",
                     "\n     no longer any \"newstyle\" argument for any",
                     "\n     of the functions.  In effect everything is now",
                     "\n     in the \"new\" style, although it is still acceptable",
                     "\n     to provide values of \"Rho\" in the \"old\" style.",
                     "\n     Such values are converted internally to the",
                     "\n     \"new\" style.\n")
	packageStartupMessage(msg)
}
