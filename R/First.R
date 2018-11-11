.onLoad <- function(lib, pkg) {
	library.dynam("hmm.discnp", pkg, lib)
}

.onAttach <- function(lib, pkg) {
	ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
        packageStartupMessage(paste(pkg, ver))
	msg <- paste("\n     This package has changed SUBSTANTIALLY from its",
                     "\n     previous release.  Read the documentation",
                     "\n     carefully.  Note in particular that the meaning of",
                     "\n     the argument \"nsim\" of the function rhmm() has",
                     "\n     changed, and a new argument \"ylengths\" now plays",
                     "\n     essentially the role previously played by \"nsim\".",
                     "\n     A new fitting method \"LM\" which uses the",
                     "\n     Levenberg-Marquardt algorithm is now available.\n")
	packageStartupMessage(msg)
}
