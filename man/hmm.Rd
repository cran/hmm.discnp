\name{hmm}
\alias{hmm}
\title{
Fit a hidden Markov model to discrete data.
}
\description{
   Effects a maximum likelihood fit of a hidden Markov model
   to discrete data where the observations come from one of a
   number of finite discrete distributions, depending on the
   (hidden) state of the Markov chain.  These distributions (the
   \dQuote{emission probabilities}) are specified non-parametrically.
   The observations may be univariate, independent bivariate,
   or dependent bivariate.  By default this function uses the EM
   algorithm.  In the univariate setting it may alternatively use
   a \dQuote{brute force} method.
}

\usage{
hmm(y, yval=NULL, par0=NULL, K=NULL, rand.start=NULL,
    method=c("EM","bf","LM","SD"), hglmethod=c("fortran","oraw","raw"),
    optimiser=c("nlm","optim"), optimMethod=NULL, stationary=cis,
    mixture=FALSE, cis=TRUE, indep=NULL, tolerance=1e-4, digits=NULL,
    verbose=FALSE, itmax=200, crit=c("PCLL","L2","Linf","ABSGRD"),
    X=NULL,keep.y=FALSE, keep.X=keep.y,
    addIntercept=TRUE, lmc=10, hessian=FALSE,...)
}
\arguments{
\item{y}{
   A vector or a list of vectors, or one or two column matrix
   (bivariate setting) or a list of such matrices; missing values
   are allowed.  If \code{y} is a vector, or list of vectors (of
   discrete data) these vectors are coerced to one column matrices.
   The entries of these vectors or matrices may be numeric or
   character and are assumed to constitute  discrete data.
} \item{yval}{
   A vector (of length \code{m}, say) of possible values for the
   data or a list of two such vectors (of lengths \code{m1} and
   \code{m2}, say, one for each of the two variates in the bivariate
   settings).  These vectors default to the sorted unique values of
   the respective variates as provided in \code{y}.  If \code{yval}
   is supplied and any value of \code{y} does not match some value
   of \code{yval}, then an error is thrown.

   The argument \code{yval} is provided so as to allow for fitting
   of models to data in which some of the data values \dQuote{of interest}
   were never observed.  The estimated emission probabilities of such
   \dQuote{never observed} values will of course be zero.
} \item{par0}{
   An optional (\emph{named}) list of starting values for the
   parameters of the model, with components \code{tpm} (transition
   probability matrix), \emph{optionally} \code{ispd} (initial state
   probability distribution) and \code{Rho}.  The object \code{Rho}
   specifies the probability that the observations take on each of
   the possible  values of the variate or variates, given the state
   of the hidden Markov chain.  See \bold{Details}.  Note that in
   the case of independent bivariate data \code{Rho} is a list
   of two matrices.  These matrices may (and in general will)
   have different row dimensions, but must have identical column
   dimensions (equal to \code{K}, the number of states; see below).

   If the model is stationary (i.e. if \code{stationary} is
   \code{TRUE}) then you should almost surely not specify the
   \code{ispd} component of \code{par0}.  If you do specify it,
   it really only makes sense to specify it to be the stationary
   distribution determined by \code{tpm} and this is a waste of
   time since this is what the code will take \code{ispd} to be if
   you leave it unspecified.

   If \code{par0} is not specified, starting values are created by
   the (undocumented) function \code{init.all()}.
}
\item{K}{
   The number of states in the hidden Markov chain; if \code{par0}
   is not specified \code{K} \emph{MUST} be; if \code{par0} is
   specified, \code{K} is ignored.

   Note that \code{K=1} is acceptable; if \code{K} is 1 then
   all observations are treated as being independent and the
   non-parametric estimate of the distribution of the observations
   is calculated in the \dQuote{obvious} way.
}
\item{rand.start}{
   Either a logical scalar or a list consisting of two logical
   scalars which must be named \code{tpm} and \code{Rho}.  If the
   former, it is converted internally into a list with entries
   named \code{tpm} and \code{Rho}, both having the same value as
   the original argument.  If \code{tpm} is TRUE then the function
   init.all() chooses entries for the starting value of \code{tpm}
   at random; likewise for \code{Rho}.  If left \code{NULL}, this
   argument defaults to \code{list(tpm=FALSE,Rho=FALSE)}.
}
\item{method}{
   Character string, either \code{"bf"}, \code{"EM"},
   \code{"LM"} or \code{"SD"} (i.e. use numerical maximisation
   via either \code{nlm()} or \code{optim()}, the EM algorithm, the
   Levenberg-Marquardt algorithm, or the method of steepest descent).
   May be abbreviated.  Currently the \code{"bf"}, \code{"LM"} and
   \code{"SD"} methods can be used only in the univariate setting,
   handle only stationary models (see below) and do not do mixtures.
}
\item{hglmethod}{
   Character string; one of \code{"fortran"}, \code{"oraw"} or
   \code{"raw"}.  May be abbreviated. This argument determines the
   procedure by which the hessian, gradient and log likelihood of
   the model and data are calculated.  If this is argument is equal
   to \code{"fortran"} (the default) then (obviously!) dynamically
   loaded fortran subroutines are used.  The other two possibilities
   effect the calculations in raw R; \code{"oraw"} (\dQuote{o}
   for \dQuote{original} uses code that is essentially a direct
   transcription of the fortran code, do-loops being replaced by
   for-loops.  With method \code{"raw"} the for-loops are eliminated
   and matrix-vector calculations are applied.  The \code{"oraw"}
   method is about 25 times slower than the \code{"fortran"} method
   and the \code{"raw"} method is (surprisingly?) even worse;
   it is more than 30 times slower.  The \dQuote{raw} methods are
   present mainly for debugging purposes and would not usually be
   used in practice.  This argument is used only if the \code{method}
   is \code{"LM"} or \code{"SD"} (and is involved only peripherally
   in the latter instance).  It is ignored otherwise.
}
\item{optimiser}{
   Character string specifying the optimiser to use when the
   \dQuote{\code{"bf"}} method of optimisation is chosen.  It should be
   one of \code{"nlm"} or \code{"optim"}, and may be abbreviated.
   Ignored unless \code{method="bf"}.
}
\item{optimMethod}{
   Character string specifying the optimisation method to be used by
   \code{\link{optim}()}.  Should be one of \code{"Nelder-Mead"},
   \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"}, or
   \code{"Brent"}.  Ignored if the \code{method} is not \code{"bf"}
   or if the optimiser is not \code{"optim"}.
}

\item{stationary}{
   Logical scalar.  If \code{TRUE} then the model is fitted under
   the stationarity assumption, i.e. that the Markov chain was in
   steady state at the time that observations commenced. In this
   case  the initial state probability distribution is estimated
   as the stationary distribution determined by the (estimated)
   transition probability matrix.  Otherwise if \code{cis} (see
   below) is \code{TRUE} the initial state probability distribution
   is estimated as the mean of the vectors of conditional
   probabilities of the states, given the observation sequences,
   at time \code{t=1}.  If \code{stationary} is \code{TRUE} and
   \code{cis} is \code{FALSE} an error is thrown.  Currently if
   the method is \code{"bf"}, \code{"LM"} or \code{"SD"}, and
   \code{stationary} is \code{FALSE}, then an error is thrown.
}
\item{mixture}{
   A logical scalar; if TRUE then a mixture model (all rows of the
   transition probability matrix are identical) is fitted rather
   than a general hidden Markov model.  Currently an error is
   thrown if \code{mixture=TRUE} and the method is
   \code{"bf"}, \code{"LM"} or \code{"SD"}.
}
\item{cis}{
   A logical scalar specifying whether there should be a
   \bold{constant} \bold{initial} \bold{state} probability
   distribution.  If \code{stationary} is \code{FALSE} and \code{cis}
   is \code{FALSE} then the initial state probability distribution
   for a given observation sequence is equal to 1 where the (first)
   maximum of the vector of conditional probabilities of the states,
   given the observation sequences, at time \code{t=1}, occurs,
   and is 0 elsewhere.  If \code{stationary} is \code{TRUE} and
   \code{cis} is \code{FALSE} an error is given.
}
\item{indep}{
   Logical scalar.  Should the bivariate model be fitted under the
   assumption that the two variables are (conditionally) independent
   give the state?  If this argument is left as \code{NULL} its
   value is inferred from the structure of \code{Rho} in \code{par0}
   if the latter is supplied.  If the data are bivariate and neither
   \code{indep} nor \code{par0} is supplied, then an error is given.
   If the data are bivariate and if the value of \code{indep}
   is inconsistent with the structure of \code{par0$Rho} then an
   error is given. If the data are univariate then \code{indep}
   is ignored.
}
\item{tolerance}{
   If the value of the quantity used for the stopping criterion
   is less than tolerance then the algorithm is considered to
   have converged.  Ignored if \code{method="bf"}.  Defaults to
   \code{1e-4}.
}
\item{digits}{
   Integer scalar.  The number of digits to which to print out
   \dQuote{progress reports} (when \code{verbose} is \code{TRUE}).
   There is a \dQuote{sensible} default (calculated from
   \code{tolerance}).  Not used if the method is \code{"bf"}.
}
\item{verbose}{
   A logical scalar determining whether to print out details of
   the progress of the algorithm.  If the method is \code{"EM"},
   \code{"LM"} or \code{"SD"} then when \code{verbose} is \code{TRUE}
   information about the convergence criteria is printed out at
   every step that the algorithm takes.  If \code{method="bf"} then
   the value of \code{verbose} determines the value of the argument
   \code{print.level} of \code{\link{nlm}()} or the value of the
   argument \code{trace} of \code{\link{optim}()}.  In the first
   case, if \code{verbose} is \code{TRUE} then \code{print.level}
   is set to 2, otherwise it is set to 0.  In the second case,
   if \code{verbose} is \code{TRUE} then \code{trace} is set to 6,
   otherwise it is set to 0.
}
\item{itmax}{
   When the method is \code{"EM"}, \code{"LM"} or \code{"SD"}
   this is the maximum number of steps that the algorithm takes.
   If the convergence criterion has not been met by the time
   \code{itmax} steps have been performed, a warning message
   is printed out, and the function stops.  A value is returned by
   the function anyway, with the logical component \code{converged} set
   to \code{FALSE}.  When \code{method="bf"} the \code{itmax} argument
   is passed to \code{nlm()} as the value of \code{iterlim}
   or to \code{optim()} as the value of \code{maxit}.  If the
   (somewhat obscure) convergence criteria of \code{nlm()} or
   \code{optim()}  have not been met by the time \code{itmax}
   \dQuote{iterations} have been performed, the algorithm ceases.
   In this case, if \code{nlm()} is used.  the value of \code{code}
   in the object returned set equal to 4 and if \code{optim()}
   is used then the value of \code{convergence} returned is set
   equal to 1.  Note that the value of \code{code}, respectively
   \code{convergence} is returned as the \code{converged} component
   of the object returned by \code{hmm()}.  A value of 1 indicates
   successful completion of the \code{nlm()} procedure.  A value of
   0 indicates successful completion of the \code{optim()} procedure.
}
\item{crit}{
   The name of the stopping criterion used.  When \code{method="EM"}
   it must be one of \code{"PCLL"} (percent change in log-likelihood;
   the default), \code{"L2"} (L-2 norm, i.e.  square root of sum of
   squares of change in coefficients), or \code{"Linf"} (L-infinity
   norm, i.e.  maximum absolute value of change in coefficients).
   When \code{method="LM"} or \code{method="SD"} there is a fourth
   possibility, namely \code{"ABSGRD"} the (maximum) absolute value
   of the gradient.  It may not be advisable to use this criterion
   in the current context (i.e. that of discrete non-parametric
   distributions).  See \bold{Warnings}.  This argument defaults
   to \code{"PCLL"}.  It is ignored if \code{method="bf"}.
   (The \code{nlm()} and \code{optim()} functions have their own
   obscure stopping criteria.)
}
\item{X}{An optional \emph{numeric} matrix, or a list of such
   matrices, of \dQuote{auxiliary} \emph{predictors}.  The use of
   such predictors is (currently, at least) applicable only in the
   univariate emissions setting.  If \code{X} is a list it must be
   of the same length as \code{y} and all entries of this list must
   have the same number of columns.  If the columns of any entry
   of the list are named, then they must be named for \emph{all}
   entries, and the column names must be the \emph{same} for all
   entries. The number of rows of each entry must be equal to the
   length of the corresponding entry of \code{y}. If \code{X} is
   a matrix then \code{y} should be a vector or one-column matrix
   (or a list with a single entry equal to such).

   There may be at most one constant column in \code{X} or the
   components thereof.  If there are \emph{any} constant columns
   there must be precisely one (in all components of \code{X}),
   it must be the first column and all of its entries must be equal
   to \code{1}.  If the columns have names, the names of this first
   column must be \code{"Intercept"}.

   Note that \code{X} (or its entries) must be a \emph{numeric}
   matrix (or must be numeric matrices) --- not data frames!  Factor
   predictors are not permitted.  It may be possible to use factor
   predictors by supplying \code{X} or its entries as the output of
   \code{\link{model.matrix}()}; this will depend on circumstances.

   The fitted coefficients that are produced when \code{X} is supplied,
   are (to put it mildly) a bit difficult to interpret.  See
   \bold{Fitted Coefficients of Auxiliary Predictors} for
   some discussion.
}
\item{keep.y}{
   Logical scalar; should the observations \code{y} be returned as
   a component of the value of this function?
}
\item{keep.X}{
   Logical scalar; should the predictors \code{X} be returned as
   a component of the value of this function?  Note that the
   value of \code{keep.X} will be silently set equal to \code{FALSE}
   unless it actually \dQuote{makes sense} to keep \code{X}.  I.e.
   unless the observations are \emph{univariate}
   and \code{X} is actually supplied, i.e. is
   not \code{NULL}.
}
\item{addIntercept}{
   Logical scalar.  Should a column of ones, corresponding to
   an intercept term, be prepended to each of the matrices in
   the list \code{X}?   If each of these matrices already has an
   initial column of ones, then setting \code{addIntercept=TRUE}
   results in an error being thrown.  If this is not the case,
   then by default an initial column of ones is added.
}

\item{lmc}{
   Numeric scalar.  The (initial) \dQuote{Levenberg-Marquardt
   correction} parameter.  Used only if \code{method="LM"},
   otherwise ignored.
}

\item{hessian}{
   Logical scalar.  Should the hessian matrix be
   returned?  This argument is relevant only if \code{method="bf"}
   (in which case it is passed along to \code{hmmNumOpt()}) and is
   ignored otherwise.  This argument should be set to \code{TRUE}
   only if you \emph{really} want the hessian matrix.  Setting it
   to \code{TRUE} causes a substantial delay between the time when
   \code{hmm()} finishes its iterations and when it actually returns
   a value.
}

\item{\dots}{
   Additional arguments passed to \code{hmmNumOpt()}.
   There is one noteworthy argument \code{useAnalGrad} which is used
   \dQuote{directly} by \code{hmmNumOpt()}.  This argument is a
   logical scalar and if it is \code{TRUE} then calls to \code{nlm()}
   or \code{optim()} are structured so that an analytic calculation
   of the gradient vector (implemented by the internal function
   \code{get.gl()} is applied.  If it is \code{FALSE} then finite
   difference methods are used to calculate the gradient vector.
   If this argument is not specified it defaults to \code{FALSE}.
   Note that the name of this argument \bold{cannot be abbreviated}.

   Other \dQuote{additional arguments} may be supplied for the
   control of \code{\link{nlm}()} and are passed on appropriately
   to \code{nlm()}.  These are used only if \code{method="bf"}
   and if \code{optimiser="nlm"}.  These \dQuote{\dots} arguments
   might typically include \code{gradtol}, \code{stepmax} and
   \code{steptol}.  They should \bold{NOT} include \code{print.level}
   or \code{iterlim}.  The former argument is automatically passed
   to \code{nlm()} as \code{0} if \code{verbose} is \code{FALSE}
   and as \code{2} if \code{verbose} is \code{TRUE}.  The latter
   argument is automatically passed to \code{nlm()} with the value
   of \code{itmax}.
}
}

\value{
   A list with components:

\item{Rho}{
   The fitted value of the data frame, list of two matrices,
   or array \code{Rho} (in the case of a univariate model, a
   bivariate independent model or a bivariate dependent model
   respectively) specifying the distributions of the observations
   (the \dQuote{emission} probabilities).
}
\item{Rho.matrix}{
   Present \emph{only} in the univariate setting.  A matrix
   whose entries are the (fitted) emission probabilities,
   row corresponding to values of the emissions and columns
   to states.  The columns sum to 1.  This component provides
   the same information as \code{Rho}, but in a more readily
   interpretable form.
}
\item{tpm}{
   The fitted value of the transition probability matrix \code{tpm}.
}
\item{stationary}{
   Logical scalar; the value of the \code{stationary} argument.
}
\item{ispd}{
   The fitted initial state probability distribution, or a matrix
   of initial state probability distributions, one (column) of
   \code{ispd} for each observation sequence.

   If \code{stationary} is \code{TRUE} then \code{ispd} is assumed
   to be the (unique) stationary distribution for the chain,
   and thereby determined by the transition probability matrix
   \code{tpm}.  If \code{stationary} is \code{FALSE} and \code{cis}
   is \code{TRUE} then \code{ispd} is estimated as the mean of the
   vectors of conditional probabilities of the states, given the
   observation sequences, at time \code{t=1}.

   If \code{cis} is \code{FALSE} then \code{ispd} is a matrix
   whose columns are the vectors of conditional probabilities of
   the states, given the observation sequences, at time \code{t=1},
   as described above.  (If there is only one observation sequence,
   then this --- one-column --- matrix is converted into a vector.)
}
\item{log.like}{
   The final (maximal, we hope!) value of the log likelihood, as
   determined by the maximisation procedure.
}
\item{grad}{
   The gradient of the log likelihood.  Present only if the
   method is \code{"LM"} or \code{"bf"} and in the latter
   case then only if the optimiser is \code{nlm()}.
}
\item{hessian}{
   The hessian of the log likelihood.  Present only if the
   method is \code{"LM"} or \code{"bf"}.
}
\item{stopCrit}{
   A vector of the (final) values of the stopping criteria, with
   names \code{"PCLL"}, \code{"L2"}, \code{"Linf"} unless the method
   is \code{"LM"} or \code{"SD"} in which case this vector has a
   fourth entry named \code{"ABSGRD"}.
}
\item{par0}{
   The starting values used by the algorithms.  Either the argument
   \code{par0}, or a similar object with either or both components
   (\code{tpm} and \code{Rho}) being created by \code{rand.start()}.
}
\item{npar}{
   The number of parameters in the fitted model.  Equal to
   \code{nispar + ntpmpar + nrhopar} where (1) \code{nispar} is
   \code{0} if \code{stationary} is \code{TRUE} and is \code{K-1}
   otherwise; (2) \code{ntpmpar} is \code{K*(K-1)} (3) \code{nrhopar}
   is
   \itemize{
   \item \code{(nrow(Rho) - K)*(ncol(Rho)-2)} for univariate models
   \item \code{K*(sum(sapply(Rho,nrow))-K)} for bivariate independent models
   \item \code{prod(dim(Rho))-K} for bivariate dependent models.
   }
}
\item{bicm}{
    Numeric scalar.  The number by which \code{npar} is multiplied
    to form the \code{BIC} criterion.  It is essentially the log
    of the number of observations.  See the code of \code{hmm()}
    for details.
}
\item{converged}{
   A logical scalar indicating whether the algorithm converged.
   If the EM, LM or steepest descent algorithm was used it simply
   indicates whether the stopping criterion was met before
   the maximum number (\code{itmax}) of steps was exceeded.
   If \code{method="bf"} then \code{converged} is based on the
   \code{code} component of the object returned by the optimiser
   when \code{\link{nlm}()} was used, or on the \code{convergence}
   component when \code{\link{optim}()} was used.  In these
   cases \code{converged} has an \emph{attribute} (\code{code}
   or \code{convergence} respectively) giving the (integer) value
   of the relevant component.

   Note that in the \code{nlm()} case a value of \code{code}
   equal to 2 indicates \dQuote{probable} convergence, and a value
   of 3 indicates \dQuote{possible} convergence.  However in this
   context \code{converged} is set equal to \code{TRUE} \emph{only}
   if \code{code} is 1.
}
\item{nstep}{
   The number of steps performed by the algorithm if the method
   was \code{"EM"}, \code{"LM"} or \code{"SD"}.  The value of
   \code{nstep} is set equal to the \code{iterations} component of
   the value returned by \code{nlm()} if \code{method="bf"}.
}

\item{prior.emsteps}{
   The number of EM steps that were taken before the method was
   switched from \code{"EM"} to \code{"bf"} or to \code{"LM"}.
   Present only in values returned under the \code{"bf"} or
   \code{"LM"} methods after a switch from \code{"EM"} and is
   equal to \code{0} if either of these methods was specified in
   the initial call (rather than arising as the result of a switch).
}
\item{ylengths}{
   Integer vector of the lengths of the observation sequences
   (number of rows if the observations are in the form of
   one or two column matrices).
}
\item{nafrac}{
   A real number between 0 and 1 or a pair (two dimensional vector)
   of such numbers.  Each number is the the fraction of missing values
   if the corresponding components of the observations.
}
\item{y}{
   An object of class \code{"tidyList"}.  It is a tidied up version
   of the observations; i.e. the observations \code{y} after the
   application of the undocumented function \code{tidyList()}.
   Present only if \code{keep.y} is \code{TRUE}.
}
\item{X}{
   An object of class \code{"tidyList"}.  It is tidied up version
   of the predictor matrix or list of predictor matrices; i.e. the
   argument \code{X} after the application of \code{tidyList()}
   (with argument \code{rp} set to \code{"predictor"}.  Present only
   if \code{X} is supplied, is an appropriate argument, and if
   \code{keep.X} is \code{TRUE}.
}
\item{parity}{
   Character string; \code{"univar"} if the data were univariate,
   \code{"bivar"} if they were bivariate.
}
\item{numeric}{
   Logical scalar; \code{TRUE} if the (original) data were numeric,
   \code{FALSE} otherwise.
}
\item{AIC}{
   The value of AIC \code{= -2*log.like + 2*npar} for the fitted
   model.
}
\item{BIC}{
   The value of BIC \code{= -2*log.like + log(nobs)*npar} for the fitted
   model.  In the forgoing \code{nobs} is the number of observations.
   This is the number of \emph{non-missing} values in \code{unlist(y)}
   in the univariate setting and one half of this number in the
   bivariate setting.
}
\item{args}{
   A list of argument values supplied.  This component is
   returned in the interest of making results reproducible.
   It is also needed to facilitate the updating of a model
   via the update method for the class \code{hmm.discnp},
   \code{\link{update.hmm.discnp}()}.

   It has components:
   \itemize{
       \item \code{method}
       \item \code{optimiser}
       \item \code{optimMethod}
       \item \code{stationary}
       \item \code{mixture}
       \item \code{cis}
       \item \code{tolerance}
       \item \code{itmax}
       \item \code{crit}
       \item \code{addIntercept}
   }
}
}

\section{Thanks}{
A massive nest of bugs was eliminated in the transition from version
3.0-8 to version 3.0-9.  These bugs arose in the context of using
auxiliary predictor variables (argument \code{X}).  The handling of
such auxiliary predictors was completely messed up.  I am grateful
to Leah Walker for pointing out the problem to me.
}

\section{Warnings}{
   The ordering of the (hidden) states can be arbitrary.  What the
   estimation procedure decides to call \dQuote{state 1} may not be
   what \emph{you} think of as being state number 1.  The ordering
   of the states will be affected by the starting values used.

   Some experiences with using the \code{"ABSGRD"} stopping
   criterion indicate that it may be problematic in the context of
   discrete non-parametric distributions.  For example a value of
   1854.955 was returned after 200 LM steps in one (non-convergent,
   of course!)  attempt at fitting a model.  The stopping criterion
   \code{"PCLL"} in this example took the \dQuote{reasonable}
   value of 0.03193748 when iterations ceased.
}

\section{Notes --- Various}{
   This function \emph{used} to have an argument \code{newstyle},
   a logical scalar (defaulting to \code{TRUE}) indicating whether
   (in the univariate setting) the emission probabilities
   should be represented in \dQuote{logistic} form.   (See
   \bold{Details}, \bold{Univariate case:}, above.)  Now the
   emission probabilities are \emph{always} represented in the
   \dQuote{logistic} form.   The component \code{Rho} of the
   starting parameter values \code{par0} may still be supplied
   as a matrix of probabilities (with columns summing to 1), but
   this component is converted (internally, silently) to the
   logistic form.

   The object returned by this function also has (in the univariate
   setting), in addition to the component \code{Rho}, a component
   \code{Rho.matrix} giving the emission probabilities in the
   more readily interpretable matrix-of-probabilities form. (See
   \bold{Value} above.)

   The package \emph{used} to require the argument \code{y} to
   be a \emph{matrix} in the case of multiple observed sequences.
   If the series were of unequal length the user was expected to
   pad them out with NAs to equalize the lengths.

   The old matrix format for multiple observation sequences was
   permitted for a while (and the matrix was internally changed into
   a list) but this is no longer allowed.  If \code{y} is indeed
   given as a matrix then this corresponds to a single observation
   sequence and it must have one (univariate setting) or two
   (bivariate setting) columns which constitute the observations
   of the respective variates.

   If \code{K=1} then \code{tpm}, \code{ispd}, \code{converged},
   and \code{nstep} are all set equal to \code{NA} in the list
   returned by this function.

   The estimate of \code{ispd} in the non-stationary setting
   is inevitably very poor, unless the number of sequences of
   observations (the length of the list \code{y}) is very large.
   We have in effect \dQuote{less than one} relevant observation for
   each such sequence.

   The returned values of \code{tpm} and \code{Rho} (or the entries
   of \code{Rho} when \code{Rho} is a list) have dimension names.
   These are formed from the argument \code{yval} if this is
   supplied, otherwise from the sorted unique values of the
   observations in \code{y}.  Likewise the returned value of
   \code{ispd} is a named vector, the names being the same as the
   row (and column) names of \code{tpm}.

   If \code{method} is equal to \code{"EM"} there may be a
   \emph{decrease} (!!!) in the log likelihood at some EM step.
   This is \dQuote{theoretically impossible} but can occur in
   practice due to an intricacy in the way that the EM algorithm
   treats \code{ispd} when \code{stationary} is \code{TRUE}.
   It turns out to be effectively impossible to maximise the expected
   log likelihood unless the term in that quantity corresponding
   to \code{ispd} is ignored (whence it \emph{is} ignored).
   Ignoring this term is \dQuote{asymptotically negligible} but
   can have the unfortunate effect of occasionally leading to a
   decrease in the log likelihood.

   If such a decrease is detected, then the algorithm terminates
   and issues a message to the effect that the decrease occurred.
   The message suggests that another method be used and that perhaps
   the results from the penultimate EM step (which are returned by
   this function) be used as starting values.

   It seems to me that it \emph{should} be the case that such a
   decrease in the log likelihood can occur only if \code{stationary}
   is \code{TRUE}.  However I have encountered instances in which
   a decrease occurred when \code{stationary} was \code{FALSE}.
   I have yet to figure out/track down what is going on here.

}
\section{Note on \code{method}}{
   If the \code{method} is \code{"EM"} it is actually possible
   for the log likelihood to \emph{decrease} at some EM step.
   This is \dQuote{impossible in an ideal world} but can happen
   to the fact the EM algorithm, as implemented in this package
   at least, cannot maximise the expected log likelihood if the
   component corresponding to the initial state probability
   distribution is taken into consideration.  This component
   should ideally be maximised subject to the constraint that
   \code{t(P)\%*\%ispd = ispd}, but this constraint seems to
   effectively impossible to impose.  Lagrangian multipliers
   don't cut it.  Hence the summand in question is ignored at
   the M-step.  This usually works alright since the summand
   is asymptotically negligible, but things can sometimes go
   wrong.  If such a decrease occurs, an error is thrown.

   In previous versions of this package, instead of throwing
   an error the \code{hmm()} function would automatically switch
   to either the \code{"bf"} or the \code{"LM"} method, depending
   whether a matrix \code{X} of auxiliary predictors is supplied,
   starting from the penultimate parameter estimates produced
   by the EM algorithm.   However this appears not to be a good
   idea; those \dQuote{penultimate estimates} appear not to be
   good starting values for the other methods.  Hence an error
   is now thrown and the user is explicitly instructed to invoke
   a different method, \dQuote{starting from scratch}.
}

\section{Fitted Coefficients of the Predictors}{
It is of course of interest to understand the meaning of the
coefficients that are fitted to the predictors in the model.
If \code{X} is supplied then the number of predictors is (as a rule)
one (for the intercept) plus the number of columns in each entry
of \code{X}.  We say \dQuote{as a rule} because, e.g., the entries
of \code{X} could each have an \dQuote{intercept} column, or the
\code{addIntercept} argument could be \code{FALSE}.  If \code{X}
is not supplied there is only one predictor, named \code{Intercept}.

The interpretation of these predictor coefficients is a bit subtle.
To get an idea of what it's all about, consider the output from
example \code{4}.  (See \bold{Examples}).  The fitted coefficients
in question are to be found in columns 3 and onward of the component
\code{Rho} of the object returned by \code{hmm()}.  In the context
of example \code{4}, this object is \code{fit.wap}.  (The suffix
\code{wap} stands for \dQuote{with auxiliary predictors}.)
\preformatted{
  fit.wap$Rho
       y state Intercept     ma.com      nh.com     bo.com
  1   lo     1 1.3810463  0.4527982 -3.27161353 -1.9563915
  2  mlo     1 0.1255631 -1.1402546 -1.37713744  0.5946980
  3    m     1 0.7356526  0.1523734 -2.70841817 -0.1794645
  4  mhi     1 0.8479798 -0.2438988 -1.12544989 -0.9650320
  5   hi     1 0.0000000  0.0000000  0.00000000  0.0000000
  6   lo     2 3.9439410 -0.8355306 -0.77702276  1.4963631
  7  mlo     2 2.6189880 -1.9373885 -0.09190623  0.8316870
  8    m     2 2.1457317 -1.7276183  0.19524655 -0.3249485
  9  mhi     2 1.8834139 -1.3760011 -0.59806309  1.2828365
  10  hi     2 0.0000000  0.0000000  0.00000000  0.0000000
}

If you multiply the matrix consisting of the predictor coefficients
(columns 3 to 6 of \code{Rho} in this instance) times a vector of
predictors you get, for each state, the \dQuote{exponential form}
of the probabilities (\dQuote{pre-probabilities}) for each of the
possible \code{y}-values, given the vector of predictors.

E.g. set \code{x <- c(1,1,0,0)}. This vector picks up the intercept
and indicates that the Malabar outfall has been commissioned,
the North Head outfall has not been commissioned, and the Bondi
Offshore outfall has not been commissioned.

Now set:
\preformatted{
    pp1 <- (as.matrix(fit.wap$Rho)[,3:6]\%*\%x)[1:5]
    pp2 <- (as.matrix(fit.wap$Rho)[,3:6]\%*\%x)[6:10]
}
Note that \code{pp1} consists of \dQuote{exponential
probabilities} corresponding to state 1, and \code{pp2}
consists of \dQuote{exponential probabilities} corresponding
to state 2.  To convert the foregoing pre-probabilities to the
actual probabilities of the \code{y}-values, we apply the ---
undocumented --- function \code{expForm2p()}:
\preformatted{
    p1 <- expForm2p(pp1)
    p2 <- expForm2p(pp2)
}

The value of \code{p1} is
\preformatted{
[1] 0.52674539 0.03051387 0.20456767 0.15400019 0.08417288
}
and that of \code{p2} is
\preformatted{
[1] 0.78428283 0.06926632 0.05322204 0.05819340 0.03503541
}
Note that \code{p1} and \code{p2} each sum to 1, as they should/must
do.  This says, e.g., that when the system is in state 2, and
Malabar has been commissioned but North Head and Bondi Offshore
have not, the (estimated) probability that \code{y} is \code{"mhi"}
(medium-high) is 0.05819340.

It may be of some interest to test the hypothesis that the predictors
have any actual predictive power at all:

\preformatted{
    fit.nap <- hmm(xxx,yval=Yval,K=2,verb=TRUE)
    # "nap" <--> no aux. preds
}

There is a bit of a problem here, in that the likelihood \emph{decreases}
at EM step 65.  (See the warning message.)

We can check on this problem by refitting using method="LM".
\preformatted{
    fit.nap.lm <- hmm(xxx,yval=Yval,par0=fit.nap,method="LM",verb=TRUE)
}
Doing so produces only a small improvement in the log likelihood
(from -1821.425 to -1820.314), so we really could have ignored the
problem. We can now do \code{anova(fit.wap,fit.nap)} which gives

\preformatted{
    $stat
    [1] 153.5491

    $df
    [1] 24

    $pvalue
    [1] 7.237102e-21
}
Thus the p-value is effectively zero, saying that in this instance
the auxiliary predictors appear to have a \dQuote{significant} impact on
the fit.

}

\details{
   \itemize{
   \item \bold{Univariate case:}
   In the univariate case the emission probabilities are specified by
   means of a data frame \code{Rho}.  The first column of \code{Rho},
   named \code{"y"}, is a factor consisting of the possible values
   of the emissions, repeated \code{K} times (where \code{K} is
   the number of states).  The second column, named \code{states},
   is a factor consisting of integer values \code{1, 2, \dots, K}.
   Each of these values is repeated \code{m} times where \code{m}
   is the length of \code{yval}.   Further columns of \code{Rho}
   are numeric and consist of coefficients of the linear predictor of
   the probabilities of the various values of \code{y}.  If \code{X}
   is \code{NULL} then \code{Rho} has only one further column named
   \code{Intercept}.

   If \code{X} is not \code{NULL} then the \code{Intercept}
   column is present only if \code{addIntercept} is \code{TRUE}.
   There as many (other, in addition to the possible \code{Intercept}
   column) numeric columns as there are columns in \code{X} or in
   the matrices in the list \code{X}.  The names of these columns
   are taken to be the column names of \code{X} or the \emph{first}
   entry of \code{X} if such column names are present.  Otherwise the
   names default to \code{V1}, \code{V2} \dots.

   The probabilities of the emissions taking on their
   various possible values are given by
   \deqn{\Pr(Y = y_i | \boldsymbol{x}, \textrm{state}=S) =
   \ell_i/\sum_{j=1}^m \ell_j}{%
   \code{Pr(Y = y_i | x, state=S) = l[i]/sum(l)}}
   where \eqn{\ell_j}{l_j} is the \eqn{j\textrm{th}}{\code{j}th}
   entry of \eqn{\boldsymbol{\beta}^{\top}\boldsymbol{x}}{%
   \code{t(beta)\%*\%x}}
   and  where in turn \eqn{\boldsymbol{x}}{\code{x}} is the vector
   of predictors and \eqn{\boldsymbol{\beta}}{\code{beta}} is the
   coefficient vector in the linear predicator that corresponds to
   \eqn{y_i} and the hidden state \eqn{S}.  For identifiability the
   vectors \eqn{\boldsymbol{\beta}}{\code{beta}} corresponding to
   the first value of \eqn{Y} (the first level of \code{Rho$y}) are
   set equal to the zero vector for all values of the state \eqn{S}.

   Note that the \code{Rho} component of the starting values
   \code{par0} may be specified as a \emph{matrix} of probabilities,
   with rows corresponding to possible values of the observations and
   columns corresponding to states.  That is the \code{Rho} component
   of \code{par0} may be provided in the form \eqn{\textrm{Rho} =
   [\rho_{ij}]}{\code{Rho = Rho[ij]}} where \eqn{\rho_{ij} = \Pr(Y = y_i
   | S = j)}{\code{Rho[i,j]= Pr(Y = y[i] | S = j)}}.  This is permissible
   as long as \code{X} is \code{NULL} and may be found to be more
   convenient and intuitive.  If the starting value for \code{Rho}
   is provided in matrix form it is (silently) converted internally
   into the data frame form, by the (undocumented) function
   \code{cnvrtRho()}.

   When argument \code{X} is \emph{not} \code{NULL}, it is
   difficult to specify a \dQuote{reasonable} value for the
   \code{Rho} component of \code{par0}.  One might try to specify
   \code{par0$Rho} in the data frame form.  The question of how
   to specify the columns of \code{par0$Rho} corresponding to the
   auxiliary predictors (columns of \code{X} or of the entries of
   \code{X}) is a thorny one.

   It is permissible in these circumstances to specify
   \code{par0$Rho} as a matrix of probabilities, just as one
   would do if \code{X} were \code{NULL}.  In this setting the
   (undocumented) function \code{checkStartVal()} converts the
   matrix of probabilities to data frame form and then appends
   columns, all of whose entries are 0, corresponding to the
   auxiliary predictors.  When \code{par0} is unspecified, the
   (undocumented) function \code{init.all()} performs similar
   construction to accommodate a non-\code{NULL} value of \code{X}.
   Whether the resulting starting value for \code{Rho} makes any
   real sense, is questionable.  However little else can be done.

   \item \bold{Independent bivariate case:} the emission
   probabilities are specified by a list of two matrices.
   In this setting \eqn{\Pr(Y_1,Y_2) = (y_{i1},y_{i2}) | S = j) =
   \rho^{(1)}_{i_1,j} \rho^{(2)}_{i_2,j}}{\code{Pr((Y[1],Y[2]) = (y[i1],y[i2])
   | S = j) = Rho.1[i1,j] * Rho.2[i2,j]}} where \eqn{R^{(k)} =
   [\rho^{(k)}_{ij}]}{\code{Rho.k}} (\eqn{k = 1,2}) are the
   two emission probability matrices.

   \item \bold{Dependent bivariate case:} the emission probabilities
   are specified by a three dimensional array.  In this setting
   \eqn{\Pr((Y_1,Y_2) = (y_{i1},y_{i2}) | S = j) = \rho_{i_1,i_2,j}}{
   \code{Pr((Y_1,Y_2) = (y[i1],y[i2]) | S = j) = Rho[i1,i2,j]}} where
   \eqn{R = [\rho_{ijk}]}{\code{Rho}} is the emission probability array.
   }

   The hard work of calculating the recursive probabilities used
   to fit the model is done by a Fortran subroutine \code{recurse}
   (actually coded in Ratfor) which is dynamically loaded.  In the
   univariate case, when \code{X} is provided, the estimation of the
   \dQuote{linear predictor} vectors \eqn{\boldsymbol{\beta}}{beta}
   is handled by the function \code{multinom()} from the \code{nnet}
   package.  Note that this is a \dQuote{Recommended} package
   and is thereby automatically available (i.e. does not have to
   be installed).
}

\references{
   Rabiner, L. R., "A tutorial on hidden Markov models and selected
   applications in speech recognition," Proc. IEEE vol. 77,
   pp. 257 -- 286, 1989.

   Zucchini, W. and Guttorp, P., "A hidden Markov model for
   space-time precipitation," Water Resources Research vol.  27,
   pp. 1917-1923, 1991.

   MacDonald, I. L., and Zucchini, W., "Hidden Markov and Other
   Models for Discrete-valued Time Series", Chapman & Hall,
   London, 1997.

   Liu, Limin, "Hidden Markov Models for Precipitation in a Region
   of Atlantic Canada", Master's Report, University of New Brunswick,
   1997.
}

\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}

\seealso{\code{\link{rhmm}()}, \code{\link{mps}()},
         \code{\link{viterbi}()}}
\examples{
# TO DO: Create one or more bivariate examples.
#
# The value of itmax in the following examples is so much
# too small as to be risible.  This is just to speed up the
# R CMD check process.
# 1.
Yval <- LETTERS[1:10]
Tpm  <- matrix(c(0.75,0.25,0.25,0.75),ncol=2,byrow=TRUE)
Rho  <- cbind(c(rep(1,5),rep(0,5)),c(rep(0,5),rep(1,5)))/5
rownames(Rho) <- Yval
set.seed(42)
xxx  <- rhmm(ylengths=rep(1000,5),nsim=1,tpm=Tpm,Rho=Rho,yval=Yval,drop=TRUE)
fit  <- hmm(xxx,par0=list(tpm=Tpm,Rho=Rho),itmax=10)
print(fit$Rho) # A data frame
print(cnvrtRho(fit$Rho)) # A matrix of probabilities
                         # whose columns sum to 1.

# 2.
# See the help for logLikHmm() for how to generate y.num.
\dontrun{
   fit.num     <- hmm(y.num,K=2,verb=TRUE,itmax=10)
   fit.num.mix <- hmm(y.num,K=2,verb=TRUE,mixture=TRUE,itmax=10)
   print(fit.num[c("tpm","Rho")])
}
# Note that states 1 and 2 get swapped.

# 3.
xxx <- with(SydColDisc,split(y,f=list(locn,depth)))
Yval <- c("lo","mlo","m","mhi","hi")
# Two states: above and below the thermocline.
fitSydCol <- hmm(xxx,yval=Yval,K=2,verb=TRUE,itmax=10)

# 4.
X <- split(SydColDisc[,c("ma.com","nh.com","bo.com")],
           f=with(SydColDisc,list(locn,depth)))
X <- lapply(X,function(x){
                 as.matrix(as.data.frame(lapply(x,as.numeric)))-1})
fit.wap <- hmm(xxx,yval=Yval,K=2,X=X,verb=TRUE,itmax=10)
# wap <--> with auxiliary predictors.

# 5.
\dontrun{ # Takes too long.
fitlm <- hmm(xxx,yval=Yval,K=2,method="LM",verb=TRUE)
fitem <- hmm(xxx,yval=Yval,K=2,verb=TRUE)
# Algorithm terminates due to a decrease in the log likelihood
# at EM step 64.
newfitlm <- hmm(xxx,yval=Yval,par0=fitem,method="LM",verb=TRUE)
# The log likelihood improves from -1900.988 to -1820.314
}

# 6.
fitLesCount <- hmm(lesionCount,K=2,itmax=10) # Two states: relapse and remission.
}
\keyword{models}
