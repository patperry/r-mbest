\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{effects}
\alias{fixef}
\alias{fixef.mhglm}
\alias{ranef}
\alias{ranef.mhglm}
\alias{VarCorr}
\alias{VarCorr.mhglm}
\alias{vcov.mhglm}
\title{
Mixed Effects
}
\description{
Get the fixed effects, random effect variances, and empirical Bayes random
effect estimates.
}
\usage{
\method{fixef}{mhglm}(object, ...)

\method{ranef}{mhglm}(object, condVar = FALSE, ...) 

\method{vcov}{mhglm}(object, ...)

\method{VarCorr}{mhglm}(x, sigma = 1, rdig = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{
    an \code{mhglm} object.
  }
  \item{sigma}{
    a factor by which to scale the random effect variance-covariance matrix.
  }
  \item{rdig}{
    ignored.
  }
  \item{condVar}{
    a logical indicating whether conditional covariance matrices for the
    random effects should be returned.
  }
  \item{\dots}{
    further arguments passed to or from other methods.
  }
}
\details{
\code{fixef} returnes the fixed effects, while \code{vcov} returns the
variance-covariance matrix of the fixed effect estimates.

\code{VarCorr} returns the random effect covariance matrix.  \code{ranef}
returns the empirical Bayes random effect estimates.

These functions behave like their counterparts in the \CRANpkg{lme4} package.
}
\seealso{
  \code{\link[lme4]{fixef}}, \code{\link[lme4]{ranef}}, \code{\link[lme4]{VarCorr}},
  from package \CRANpkg{lme4}.
}
\keyword{models}
