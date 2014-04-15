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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under \sQuote{Details}.}

  \item{family}{a description of the error distribution and link
    function to be used in the model. This can be a character string
    naming a family function, a family function or the result of a call
    to a family function.  (See \code{\link{family}} for details of
    family functions.)}

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{hglm.fast} is called.}

  \item{weights}{an optional vector of \sQuote{prior weights} to be used
    in the fitting process.  Should be \code{NULL} or a numeric vector.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}

  \item{start}{starting values for the fixed effect parameters in the linear predictor.}

  \item{etastart}{starting values for the linear predictor.}

  \item{mustart}{starting values for the vector of means.}

  \item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.
    This should be \code{NULL} or a numeric vector of length equal to
    the number of cases.  One or more \code{\link{offset}} terms can be
    included in the formula instead or as well, and if more than one is
    specified their sum is used.  See \code{\link{model.offset}}.}

  \item{control}{a list of parameters for controlling the fitting
    process.  For \code{hglm.fast.fit} this is passed to
    \code{\link{hglm.fast.control}}.}

  \item{model}{a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value.}

  \item{method}{the method to be used in fitting the model.  The default
    method \code{"hglm.fast.fit"} uses moment-based estimates;
    the alternative \code{"model.frame"} returns the model frame
    and does no fitting.

    User-supplied fitting functions can be supplied either as a function
    or a character string naming a function, with a function which takes
    the same arguments as \code{hglm.fast.fit}.  If specified as a character
    string it is looked up from within the \pkg{mbest} namespace.
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

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{intercept}{logical. Should an intercept be included in the
    \emph{null} model?}
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
  \code{reterms | g1 /.../gQ} are not supported.
}
\value{
  \code{hglm.fast} returns an object of class inheriting from \code{"hglm"}.
  See later in this section.  If a non-standard \code{method} is used, the
  object will also inherit from the class (if any) returned by that function.

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
  If the moment-based random effect covariance is not positive-definite, then
  a warning will be issued, and a projection of the estimate to the
  positive-semidefinite cone will be used instead.
}

\examples{
}
\keyword{models}
\keyword{regression}
