\name{lesionCount}
\alias{lesionCount}
\docType{data}
\title{Multiple sclerosis lesion counts for three patients.}
\description{
   Lesion counts for three multiple sclerosis patients.
   The counts were obtained by magnetic resonance imaging,
   and were observed at monthly intervals.
}
\usage{
  lesionCount
}
\format{
    A list with three components each component being the sequence
    of counts for a given patient and consisting of a vector with
    non-negative integer entries.
}

\source{
   The data were originally studied by Albert et al., (1994).
   They are were also analyzed by Altman and Petkau (2005).
   The data were kindly provided by Prof. Altman.
}
\references{
Albert, P. S., McFarland, H. F., Smith, M. E., and Frank, J. A.
Time series for modelling counts from a relapsing-remitting disease:
application to modelling disease activity in multiple sclerosis.
\emph{Statistics in Medicine} \bold{13} (1994) 453--466.

Altman, Rachel MacKay, and Petkau, A. John.  Application
of hidden Markov models to multiple sclerosis lesion count data.
\emph{Statistics in Medicine} \bold{24} (2005) 2335--2344.
}

\section{Modelling}{
The hidden Markov models applied to these data by Albert et al.
and by MacKay and Petkau were much more complex and elaborate than
those fitted in the examples in this package.  See the references
for details.
}

\keyword{datasets}
