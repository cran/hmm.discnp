.onAttach <- function(lib, pkg) {
	ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
        packageStartupMessage(paste(pkg, ver))
	msg <- paste("\n     A massive nest of bugs, that were present in the",
                     "\n     previous release of this package, has been",
                     "\n     eliminated in the current release.  See the help",
                     "\n     for hmm() to get a bit more detail.  Thanks to",
                     "\n     Leah Walker for pointing out the problem.\n")
	packageStartupMessage(msg)
}
