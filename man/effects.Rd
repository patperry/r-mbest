\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{effects}
\alias{fixef}
\alias{fixef.mhglm}
\alias{fixef.mhglm_ml}
\alias{ranef}
\alias{ranef.mhglm}
\alias{ranef.mhglm_ml}
\alias{VarCorr}
\alias{VarCorr.mhglm}
\alias{VarCorr.mhglm_ml}
\alias{vcov.mhglm}
\alias{vcov.mhglm_ml}
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

\method{VarCorr}{mhglm}(x, sigma = 1, ...)

\method{fixef}{mhglm_ml}(object, ...)

\method{ranef}{mhglm_ml}(object, condVar = FALSE, ...) 

\method{vcov}{mhglm_ml}(object, ...)

\method{VarCorr}{mhglm_ml}(x, sigma = 1, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{
    an \code{mhglm} object.
  }
  \item{sigma}{
    a factor by which to scale the random effect variance-covariance matrix.
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

These functions behave like their counterparts in the \CRANpkg{nlme} package.
}
\seealso{
  \code{\link[nlme]{fixef}}, \code{\link[nlme]{ranef}}, \code{\link[nlme]{VarCorr}},
  from package \CRANpkg{nlme}.
}
\keyword{models}
