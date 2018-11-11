\name{hydroDat}
\alias{linLandFlows}
\alias{ftLiardFlows}
\alias{portMannFlows}
\alias{portMannSedLoads}
\alias{portMannSedCon}
\docType{data}
\title{
  Canadian hydrological data sets.
}
\description{
  Five data sets obtained from the \dQuote{HYDAT} database,
  Environment and Climate Change Canada's database of historical
  hydrometric data.  The data were obtained using the \code{tidyhydat}
  package.  The data have been trimmed so that there are no gaps in
  the observation dates and are presented in \dQuote{raw} form and in
  discretised form as deciles of the residuals (difference between
  raw values and the daily mean over years).
}
\usage{
    data("linLandFlows")
    data("ftLiardFlows")
    data("portMannFlows")
    data("portMannSedLoads")
    data("portMannSedCon")
}

\format{
  Data frames with observations on the following 3 variables.
  \describe{
    \item{\code{Date}}{Dates on which observations were made.}
    \item{\code{Value}}{Numeric vector of observation values.}
    \item{\code{mean}}{The mean over years of \code{Value}.}
    \item{\code{resid}}{The difference \code{Value - mean}.}
    \item{\code{deciles}}{A factor with levels \code{d1}, \dots, \code{d10},
     which are the deciles of the variable \code{resid}}
  }
}
\details{
   The variable \code{mean} was calculated as follows:
   \preformatted{
        yday <- as.POSIXlt(X$Date)$yday
        mn   <- tapply(X$Value,yday,mean,na.rm=TRUE)
        mean <- mn[as.character(yday)]
   }
   where \code{X} is the data set being processed.

   The data set \code{linLandFlows} originally consisted of 2008 observations;
   there were 1980 observations after \dQuote{trimming}.
   The data set \code{ftLiardFlows} originally consisted of 22364 observations;
   there were 11932 observations after \dQuote{trimming}.
   The data set \code{portMannFlows} originally consisted of 6455 observations;
   there were 3653 observations after \dQuote{trimming}.
   The data set \code{portMannSedLoads} consists of 2771 observations;
   no observations were trimmed.
   The data set \code{portMannSedCon} consists of 4597 observations;
   no observations were trimmed.

   The units of the \dQuote{Flows} variables are cubic metres per
   second (\eqn{m^3/s}); the units of \dQuote{portMannSedLoads}
   are tonnes; the units of \dQuote{portMannSedCon} are milligrams
   per litre (mg/l).

   The \dQuote{linLandFlows} data were obtained at the Lindberg
   Landing hydrometric station on the Liard River in the Northwest
   Territories of Canada.  The \dQuote{ftLiardFlows} data were
   obtained at the Fort Liard hydrometric station on the Liard River
   in the Northwest Territories of Canada.  The \dQuote{portMann}
   data were obtained at the hydrometric station located at the
   Port Mann pumping station on the Fraser River in the Province
   of British Columbia in Canada.
}
\source{
  Environment and Climate Change Canada's database \dQuote{HYDAT},
  a database of historical hydrometric data.  The data were obtained
  vis the \code{tidyhydat} package, which is available from \dQuote{CRAN},
  \code{https://cran.r-project.org}
}
\examples{
fit <- hmm(linLandFlows$deciles,K=4,itmax=10)
}
\keyword{datasets}
