\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{mhglm}
\alias{mhglm}
\alias{mhglm.fit}
\title{
Fitting Moment Hierarchical Generalized Linear Models
}
\description{
\code{mhglm} is used to fit a moment hierarchical generalized linear model.
}
\usage{
mhglm(formula, family = gaussian, data, weights, subset,
      na.action, start = NULL, etastart, mustart, offset,
      control = list(), model = TRUE, method = "mhglm.fit",
      x = FALSE, z = FALSE, y = TRUE, group = TRUE,
      contrasts = NULL)

mhglm.fit(x, z, y, group, weights = rep(1, nobs),
          start = NULL, etastart = NULL, mustart = NULL,
          offset = rep(0, nobs), family = gaussian(), 
          control = list(), intercept = TRUE)
}
\arguments{
  \item{formula, family, data, weights, subset, na.action, start, etastart,
    mustart, offset, model, contrasts, intercept}{These arguments
    are analogous to the similarly-named arguments for the \code{\link{glm}} and
    \code{\link{glm.fit}} functions.
  }

  \item{control}{a list of parameters for controlling the fitting
    process.  For \code{mhglm.fit} this is passed to
    \code{\link{mhglm.control}}.}

  \item{method}{the method to be used in fitting the model.  The default
    method \code{"mhglm.fit"} uses moment-based estimates;
    the alternative \code{"model.frame"} returns the model frame
    and does no fitting.
  }

  \item{x, z, y, group}{For \code{mhglm}:
    logical values indicating whether the response vector, model matrices,
    and grouping factor used in the fitting process should be returned as
    components of the returned value.

    For \code{mhglm.fit}: \code{x} is a fixed effect design matrix of
    dimension \code{n * p}, \code{z} is a random effect design matrix of
    dimension \code{n * q}, \code{y} is a vector of observations of length
    \code{n}, and \code{group} is a grouping factor vector of length
    \code{n}.
  }
}
\details{
  These functions are analogues of \code{\link{glm}} and
  \code{\link{glm.fit}}, meant to be used for fitting hierarchical
  generalized linear models.  A typical predictor has the form
  \code{response ~ terms + (reterms | group)} where
  \code{response} is the (numeric) response vector, \code{terms} is a
  series of terms which specifies a linear predictor for
  \code{response}, \code{reterms} is a series of terms with random
  coefficients (effects), and \code{group} is a grouping factor; observations
  with the same grouping factor share the same random effects.

  Currently, only one random effect term is allowed, along with a single
  level of hierarchy; random effect terms of the form
  \code{reterms | g1/.../gQ} are not supported.
}
\value{
  \code{mhglm} returns an object of class inheriting from \code{"mhglm"}.

  The function \code{\link{summary}} can be used to obtain or print a summary
  of the results.

  The generic accessor functions \code{\link{fixef}}, \code{\link{ranef}},
  \code{\link{VarCorr}}, \code{\link{sigma}}, \code{fitted.values} and
  \code{residuals} can be used to extract various useful features of the
  value returned by \code{mhglm}.
}
\seealso{
  \code{\link{terms.mhglm}}, \code{\link{model.matrix.mhglm}}, and
  \code{\link{predict.mhglm}} for \code{mhglm} methods, and the
  generic functions \code{\link{fitted.values}}, \code{\link{residuals}},
  \code{\link{summary}}, \code{\link{vcov}}, and \code{\link{weights}}.

  Generic functions \code{\link{fixef}}, \code{\link{ranef}},
  \code{\link{VarCorr}}, and \code{\link{sigma}} for features
  related to mixed effect models.

  \code{\link{glmer}} (package \CRANpkg{lme4}) for
  fitting generalized linear mixed models with likelihood-based estimates.
}
\references{
  Perry, P. O. (2015)
  "Fast Moment-Based Estimation for Hierarchical Models",
  Preprint.
}
\author{
  Patrick O. Perry
}
\note{
  If the moment-based random effect covariance is not positive-semidefinite, then
  a warning will be issued, and a projection of the estimate to the
  positive-semidefinite cone will be used instead.
}

\examples{
library(lme4)
## The following examples are adapted from lme4:
(fm1 <- mhglm(Reaction ~ Days + (Days | Subject), gaussian, sleepstudy))

(gm <- mhglm(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial))
}
\keyword{models}
\keyword{regression}
