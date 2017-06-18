\name{mhglm.control}
\alias{mhglm.control}
\alias{mhglm_ml.control}
\title{
Auxiliary for Controlling Moment Heirarchical GLM Fitting
}
\description{
Auxiliary function for \code{\link{mhglm}} fitting. Typically only used
internally by \code{\link{mhglm.fit}}, but may be used to construct a
control argument to either function.
}
\usage{
mhglm.control(standardize = TRUE, steps = 1, parallel = FALSE, diagcov = FALSE,
              fit.method = "firthglm.fit", fit.control = list(...), ...)

mhglm_ml.control(standardize = FALSE, steps = 1, parallel = FALSE, diagcov = FALSE,
              fit.method = "firthglm.fit", fit.control = list(...), ...)

}
\arguments{
  \item{standardize}{
    logitcal indicating if predictors should be standardized before
    moment-based fitted
}
  \item{steps}{
    number of refinement steps
}
  \item{parallel}{
    fit the group-specific estimates in parallel rather than sequentially
}
  \item{diagcov}{
    estimate random effect covairance matrix with diagonal approximation 
}
  \item{fit.method}{
    method for obtaining group-specific effect estimates
}
  \item{fit.control}{
    control parameters for \code{fit.method}
}
  \item{\dots}{
    arguments to be used to form the \code{fit.control} argument if
    it is not supplied directly.
}
}
\details{
Setting \code{standardize = TRUE} ensures that the procedure is equivariant,
and generally leads to better estimation performance.
Right now \code{standardize = TRUE} is not allowed for \code{mhglm_ml}.

The \code{steps} argument gives the number of refinement steps for the moment
based parameters.  In each step, the previous fixed effect and random effect
covariance matrix estimates are used to weight the subpopulation-specific
effect estimates.  In principle, higher values of \code{steps} could lead to
more accurate estimates, but in simulations, the differences are negligible.
}
\value{
A list with components named as the arguments.
}
\seealso{
\code{\link{mhglm.fit}}, the fitting procedure used by
\code{\link{mhglm}}.

\code{\link{firthglm.fit}}, the default subpopulation-specific fitting method.
}
\examples{
library(lme4) # for cbpp data

# The default fitting method uses Firth's bias-corrected estimates
(gm.firth <- mhglm(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   data = cbpp, family = binomial,
                   control=mhglm.control(fit.method="firthglm.fit")))

# Using maximum likelihood estimates is less reliable
(gm.ml <- mhglm(cbind(incidence, size - incidence) ~ period + (1 | herd),
                data = cbpp, family = binomial,
                control=mhglm.control(fit.method="glm.fit")))
}
\keyword{optimize}
\keyword{models}
