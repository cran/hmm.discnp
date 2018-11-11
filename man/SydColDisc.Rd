\name{SydColDisc}
\alias{SydColDisc}
\docType{data}
\title{Discretised version of coliform counts in sea-water samples}
\description{
   Discretised version of counts of faecal coliform bacteria in sea water
   samples collected at seven locations near Sydney NSW, Australia.
   There were four \dQuote{controls}: Longreef, Bondi East, Port
   Hacking \dQuote{50}, and Port Hacking \dQuote{100} and three
   \dQuote{outfalls}: Bondi Offshore, Malabar Offshore and North
   Head Offshore.  At each location measurements were made at four
   depths: 0, 20, 40, and 60 meters.  A large fraction of the counts
   are missing values.
}
\usage{SydColDisc}
\format{
  A data frame with 5432 observations on the following 6 variables.
  \describe{
    \item{\code{y}}{A factor consisting of a discretisation of counts
    of faecal coliform count bacteria in sea water samples.  The
    original measures were obtained by a repeated dilution process.}
    The data were discretised using the \code{\link{cut}()}
    function with breaks given by \code{c(0,1,5,25,200,Inf)} and
    labels equal to \code{c("lo","mlo","m","mhi","hi")}.
    \item{\code{locn}}{a factor with levels \code{Longreef},
    \code{Bondi East}, \code{Port Hacking 50}, \code{Port Hacking
    100}, \code{Bondi Offshore}, \code{Malabar Offshore} and
    \code{North Head Offshore}.}
    \item{\code{depth}}{a factor with levels \code{0} (0 metres),
    \code{20} (20 metres), \code{40} (40 metres), \code{60} (60 metres).}
    \item{\code{ma.com}}{A factor with levels \code{no} and \code{yes},
    indicating whether the Malabar sewage outfall had been commissioned.}
    \item{\code{nh.com}}{A factor with levels \code{no} and \code{yes},
    indicating whether the North Head sewage outfall had been commissioned.}
    \item{\code{bo.com}}{A factor with levels \code{no} and \code{yes},
    indicating whether the Bondi Offshore sewage outfall had been commissioned.}
  }
}

\details{
   The observations corresponding to each location-depth combination
   constitute a (discrete valued) time series.  The sampling
   interval is ostensibly 1 week; distinct time series are ostensibly
   synchronous.  The measurements were made over a 194 week period.
   Due to exigencies of weather, the unreliabitity of boats and
   other factors the collection times were actually highly irregular
   and have been rounded to the neares week.  Often no sample was
   obtained at a given site within a week of the putative collection
   time, in which the observed count is given as a missing value. In
   fact over \bold{75\%} of the counts are missing.  See Turner et
   al. (1998) for more detail.
}

\source{
  The original data were kindly supplied by Geoff Coade, of the
  New South Wales Environment Protection Authority (Australia)
}

\references{
  T. Rolf Turner, Murray A. Cameron, and Peter J. Thomson.  Hidden
  Markov chains in generalized linear models.  \emph{Canadian J. Statist.}
  \bold{26} (1998) 107 -- 125.
  
  Rolf Turner.  Direct maximization of the likelihood of a hidden
  Markov model.  \emph{Comp. Statist. Data Anal.} \bold{52} (2008)
  4147--4160.
}

\section{Modelling}{
The hidden Markov models applied in the paper Turner et
al. (1998) and in the paper Turner (2008) used a numeric version
of the response in this data set.  The numeric response was
essentially a square root transformation of the original data,
and the resulting values were modelled in terms of a Poisson
distribution.  See the references for details.
}

\keyword{datasets}
