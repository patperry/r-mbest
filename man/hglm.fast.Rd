\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{hglm.fast}
\alias{hglm.fast}
\alias{hglm.fast.fit}
\title{
Fitting Hierarchical Generalized Linear Models
}
\description{
\code{hglm.fast} is used to fit a hierarchical generalized linear model with
fast moment-based estimates.
}
\usage{
hglm.fast(formula, family = gaussian, data, weights, subset,
          na.action, start = NULL, etastart, mustart, offset,
          control = list(), model = TRUE, method = "hglm.fast.fit",
          x = FALSE, z = FALSE, y = TRUE, group = TRUE,
          contrasts = NULL)

hglm.fast.fit(x, z, y, group, weights = rep(1, nobs),
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
    process.  For \code{hglm.fast.fit} this is passed to
    \code{\link{hglm.fast.control}}.}

  \item{method}{the method to be used in fitting the model.  The default
    method \code{"hglm.fast.fit"} uses moment-based estimates;
    the alternative \code{"model.frame"} returns the model frame
    and does no fitting.
  }

  \item{x, z, y, group}{For \code{hglm.fast}:
    logical values indicating whether the response vector, model matrices,
    and grouping factor used in the fitting process should be returned as
    components of the returned value.

    For \code{hglm.fast.fit}: \code{x} is a fixed effect design matrix of
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
  \code{hglm.fast} returns an object of class inheriting from \code{"hglm"}.

  The function \code{\link{summary}} (i.e., \code{\link{summary.hglm}}) can
  be used to obtain or print a summary of the results.

  The generic accessor functions \code{\link[nlme]{fixef}}, \code{\link[nlme]{ranef}},
  \code{fitted.values} and \code{residuals} can be used to
  extract various useful features of the value returned by \code{hglm.fast}.
}
\seealso{
  \code{\link{terms.hglm}}, \code{\link{model.matrix.hglm}}, and
  \code{\link{predict.hglm}} for \code{hglm} methods, and the
  generic functions \code{\link{residuals}},
  \code{\link{summary}}, \code{\link{vcov}}, \code{\link{weights}}.

  Generic functions \code{\link[nlme]{fixef}}, \code{\link[nlme]{ranef}},
  \code{\link[nlme]{VarCorr}}, and \code{\link[lme4]{sigma}} for features
  related to mixed effect models.

  \code{\link[lme4]{glmer}} (package \CRANpkg{lme4}) for
  fitting generalized linear mixed models with likelihood-based estimates.
}
\references{
  Perry, P. O. (2014)
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
}
\keyword{models}
\keyword{regression}
