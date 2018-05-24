\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{mhglm}
\alias{mhglm}
\alias{mhglm.fit}
\alias{mhglm_ml}
\alias{mhglm_ml.fit}
\title{
Fitting Moment Hierarchical Generalized Linear Models
}
\description{
\code{mhglm} is used to fit a moment hierarchical generalized linear model of one level.
\code{mhglm_ml} is used to fit a moment hierarchical generalized linear model of arbitrary number of levels 
(including one level).
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
          control = list(), intercept = TRUE, dispersion = NULL)

mhglm_ml(formula, family = gaussian, data, weights, subset,
         na.action, start = NULL, etastart, mustart, offset,
         control = list(), model = TRUE, method = "mhglm_ml.fit",
         x = FALSE, z = FALSE, y = TRUE, group = TRUE,
         contrasts = NULL)

mhglm_ml.fit(x, z, y, group, weights = rep(1, nobs),
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

  \item{x, z, y, group}{ For \code{mhglm} and \code{mhglm_ml}:
      logical values indicating whether the response vector, model matrices,
      and grouping factor used in the fitting process should be returned as
      components of the returned value.

    For \code{mhglm.fit}: \code{x} is a fixed effect design matrix of
    dimension \code{n * p}, \code{z} is a random effect design matrix of
    dimension \code{n * q}, \code{y} is a vector of observations of length
    \code{n}, and \code{group} is a grouping factor vector of length
    \code{n}.

    For \code{mhglm_ml.fit}: \code{x} is a fixed effect design matrix of
    dimension \code{n * p}, \code{z} is a list of L elements, with L
    the depth of nested hierarchies, each element of \code{z} 
    is a random effect design matrix of dimension \code{n * q_i}, with
    q_i the feature dimension on tree depth i, 
    \code{y} is a vector of observations of length \code{n}, 
    and \code{group} is a list of L elements (same L as \code{z}), 
    each element of \code{group} is a grouping factor vector of length \code{n}.

  }

  \item{dispersion}{If \code{NULL}, will estimate from data; otherwise use this argument as 
      dispersion parameter.}


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

  \code{mhglm} and \code{mhglm.fit} only allow one random effect term,
  along with a single level of hierarchy. 
  \code{mghlm_ml} and \code{mhglm_ml.fit} allow multiple random effect terms so 
  long as levels of random effects are hierarchically nested. If the random 
  effect design matrices are the same for each level, a predictor has the form
  \code{response ~ terms + (reterms | g1/.../gQ)}. If the random effects design
  matrices differ from level to level, colons are used to delineate the nesting
  structure; for example, 
  \code{response ~ fe + (re1 | g1) + (re2 | g2:g1) + (re3 | g3:g2:g1)}.

  \code{mhglm} allows \code{||} in the formula
  \code{response ~ terms + (reterms || group)} to indicate that random effects
  are independent, that is the random effects covariance matrix has non-zero 
  value only on its diagonal.
  \code{mhglm_ml} currently does not support \code{||}, to indicate indpenedent
  random effects, set \code{control=list(diagcov = TRUE)}.
}

\value{
  \code{mhglm} returns an object of class inheriting from \code{"mhglm"}.
  \code{mhglm_ml} returns an object of class inheriting from \code{"mhglm_ml"}.

  The function \code{\link{summary}} can be used to obtain or print a summary
  of the results.

  The generic accessor functions \code{\link{fixef}}, \code{\link{ranef}},
  \code{\link{VarCorr}}, \code{\link{sigma}}, \code{fitted.values} and
  \code{residuals} can be used to extract various useful features of the
  value returned by \code{mhglm} or \code{mhglm_ml}.
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
P. O. Perry (2017)
"Fast moment-based estimation for hierarchical models."

N. Zhang, K. Schmaus, and P. O. Perry (2018)
"Fitting deeply-nested hierarchical models to a large book review dataset using moment-based estimators."
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

(fm2 <- mhglm(Reaction ~ Days + (Days || Subject), gaussian, sleepstudy))

(gm <- mhglm(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial))

## The following examples are for multilevel models
g_data <- mhglm_sim(n = 30, m_per_level = c(10, 5, 2), sd_intercept = c(1, 1, 1),
                    sd_slope = c(1, 1, 1), family = "gaussian", seed = 12345)

(m1 <- mhglm_ml(y ~ 1 + x + (1 + x | g1/g2/g3), gaussian, g_data))
# or equivalent
(m2 <- mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g2:g1) + (1 + x | g3:g2:g1),
                gaussian, g_data))
}
\keyword{models}
\keyword{regression}
