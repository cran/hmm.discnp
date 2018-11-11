\name{ccprSim}
\alias{ccprSim}
\docType{data}
\title{
Simulated monocyte counts and psychosis symptoms.
}
\description{
Discretised values of monocyte counts, and ratings of level of
psychosis simulated from a model fitted to a data set consisting of
observations made on a number of patients from the Northern District
Health Board system.  The real data must be kept confidential due
to ethics constraints.
}
\usage{
    data("ccprSim")
}
\format{
    The object \code{ccprSim} is a list of length 1258.  Each entry
    of this list is to be considered to correspond to an individual
    subject.  The entries consist of matrices having two columns
    named \code{cellCount} and \code{psychosisRating}.  The number
    of rows of these matrices varies from entry to entry of the list
    (i.e. from subject to subject).

    Most of the entries of these matrices are \code{NA}.  The entries
    are temporally ordered and correspond to the number of weeks
    from the start of observation.  Observations in the real
    data set were made only when the patient in question visted a
    physician and so weeks in which no visit was made resulted in an
    \dQuote{observation} of \code{NA}.  The object \code{ccprSim}
    was simulated in such a way as to imitate this characteristic.
    The fraction of missing observations in each variate (i.e.
    \code{cellCount} and \code{psychosisRating} is roughly
    commensurate with the corresponding fractions in the real data.

    The values in the first column of each matrix (the
    \code{cellCount} column consist of integers from 1 to 5 and
    are to be interpreted as indicators of cell counts in units of
    \eqn{10^9}{1e9} cells per litre, discretised according to the
    following scale:
    \itemize{
    \item \eqn{0.0 \leq c \leq 0.3}{0.0 <= c <= 0.3}
          \eqn{\leftrightarrow}{\code{<-->}} 1
    \item \eqn{0.3 < c \leq 0.5}{0.3 < c <= 0.5}
          \eqn{\leftrightarrow}{\code{<-->}} 2
    \item \eqn{0.5 < c \leq 0.7}{0.5 < c <= 0.7}
          \eqn{\leftrightarrow}{\code{<-->}} 3
    \item \eqn{0.7 < c \leq 1.0}{0.7 < c <= 1.0}
          \eqn{\leftrightarrow}{\code{<-->}} 4
    \item \eqn{1.0 < c \leq 2.0}{1.0 < c <= 2.0}
          \eqn{\leftrightarrow}{\code{<-->}} 5
     }
     where \eqn{c} represents \dQuote{count}.

    The values in the second column of each matrix (the
    \code{psychosisRating} column consist of integers from 0 to
    4 and are to be interpreted as indicators of a physician's
    assessment of the level of pschosis of the patient.  A value of
    0 corresponds to \dQuote{no symptoms}; a value of 4 corresponds
    to \dQuote{severe}.

    The question of essential interest in respect of the real data
    was \dQuote{Is there any association between the cell count
    values and the psychosis ratings?}  More specifically it was
    \dQuote{Can the level of psychosis be \emph{predicted} from
    the cell counts?}
}

\source{
  The real data, on the basis of which these data were simulated, were
  supplied by Dr. Jonathan Williams, Northern District Health Board.
}
\examples{
\dontrun{
# Takes a long time:
    fit <- hmm(ccprSim,K=2,indep=FALSE,itmax=500,verbose=TRUE)
# Throws an error; inadequate data:
    subsetFit <- hmm(ccprsim[1021:1025],K=2,indep=FALSE,itmax=5)
}
}
\keyword{datasets}
