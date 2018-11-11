\encoding{UTF-8}
\name{weissData}
\alias{weissData}
\alias{Bovine}
\alias{Cryptosporidiosis}
\alias{Downloads}
\alias{EricssonB_Jul2}
\alias{FattyLiver}
\alias{FattyLiver2}
\alias{goldparticle380}
\alias{Hanta}
\alias{InfantEEGsleepstates}
\alias{IPs}
\alias{LegionnairesDisease}
\alias{OffshoreRigcountsAlaska}
\alias{PriceStability}
\alias{Strikes}
\alias{WoodPeweeSong}
\docType{data}
\title{
   Data from \dQuote{An Introduction to Discrete-Valued Time Series}
}
\description{
   Data sets from the book
   \dQuote{An Introduction to Discrete-Valued Time Series} by
   Christian H. \enc{Weiß}{Weiss}.
}
%% To get the ß symbol use <compose key> s s.
\usage{
    data(Bovine)
    data(Cryptosporidiosis)
    data(Downloads)
    data(EricssonB_Jul2)
    data(FattyLiver)
    data(FattyLiver2)
    data(goldparticle380)
    data(Hanta)
    data(InfantEEGsleepstates)
    data(IPs)
    data(LegionnairesDisease)
    data(OffshoreRigcountsAlaska)
    data(PriceStability)
    data(Strikes)
    data(WoodPeweeSong)
}
\format{
  \itemize{
  \item \code{Bovine}
  A character vector of length 8419.
  \item \code{Cryptosporidiosis}
  A numeric (integer) vector of length 365.
  \item \code{Downloads}
  A numeric (integer) vector of length 267.
  \item \code{EricssonB_Jul2}
  A numeric (integer) vector of length 460.
  \item \code{FattyLiver2}
  A numeric (integer) vector of length 449.
  \item \code{FattyLiver}
  A numeric (integer) vector of length 928.
  \item \code{goldparticle380}
  A numeric (integer) vector of length 380.
  \item \code{Hanta}
  A numeric (integer) vector of length 52.
  \item \code{InfantEEGsleepstates}
  A character vector of length 107.
  \item \code{IPs}
  A numeric (integer) vector of length 241.
  \item \code{LegionnairesDisease}
  A numeric (integer) vector of length 365.
  \item \code{OffshoreRigcountsAlaska}
  A numeric (integer) vector of length 417.
  \item \code{PriceStability}
  A numeric (integer) vector of length 152.
  \item \code{Strikes}
  A numeric (integer) vector of length 108.
  \item \code{WoodPeweeSong}
  A numeric (integer) vector of length 1327.
  }
}

\details{
   For detailed information about each of these data sets, see the
   book cited in the \bold{References}.

   Note that the data sets \code{Cryptosporidiosis}
   and \code{LegionnairesDisease} are actually
   called \cr \code{Cryptosporidiosis_02-08} and
   \code{LegionnairesDisease_02-08} in the given reference.
   The \cr \dQuote{suffixes} were removed since the minus sign causes
   problems in a variable name in \code{R}.
}

\source{
These data sets were kindly provided by Prof. Christian
H. \enc{Weiß}{Weiss}.  The package author is also pleased
to acknowledge the kind permission granted by Prof. Kurt
\enc{Brännäs}{Brannas} (Professor Emeritus of Economics at
\enc{Umeå}{Umea} University) to include the Ericsson time series
data set (\code{EricssonB_Jul2}).
}
%% Email: weissc@hsu-hh.de
\references{Christian H. \enc{Weiß}{Weiss} (2018).  \emph{An
Introduction to Discrete-Valued Time Series}.  Chichester: John
Wiley & Sons.
}
\examples{
\dontrun{
fit1 <- hmm(WoodPeweeSong,K=2,verbose=TRUE)
# EM converges in 6 steps --- suspicious.
set.seed(321)
fit2 <- hmm(WoodPeweeSong,K=2,verbose=TRUE,rand.start=list(tpm=TRUE,Rho=TRUE))
# 52 steps --- note the huge difference between fit1$log.like and fit2$log.like!
set.seed(321)
fit3 <- hmm(WoodPeweeSong,K=2,verbose=TRUE,method="bf",
            rand.start=list(tpm=TRUE,Rho=TRUE))
# log likelihood essentially the same as for fit2
}
}
\keyword{datasets}
