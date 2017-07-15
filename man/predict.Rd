\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{predict}
\alias{predict.mhglm}
\alias{predict.mhglm_ml}
\alias{sigma}
\alias{sigma.mhglm}
\alias{sigma.mhglm_ml}
\title{
Prediction
}
\description{
\code{predict} gives empirical Bayes predictions of the response, while
\code{sigma} gives the dispersion parameter.
}
\usage{
\method{predict}{mhglm}(object, newdata = NULL, type = c("link", "response"),
      se.fit = FALSE, na.action = na.pass, ...)

\method{sigma}{mhglm}(object, ...)

}
\arguments{
  \item{object}{
    an \code{mhglm} object
  }
  \item{newdata, type, se.fit, na.action}{
    these arguments behave as in \code{\link{predict.glm}}.  See Details, below.
  }
  \item{\dots}{
    further arguments passed to or from other methods.
  }
}
\details{
The \code{predict} function gives empirical Bayes posterior mean estimates of
response values.  If \code{se.fit = TRUE}, then the conditional variances of
the random effects are used along with the fixed effect variance-covariance
matrix to estimate the standard errors.

The \code{sigma} function gives the square root of the dispersion parameter
or the model; for linear models, this is the error standard deviation.
}
\seealso{
\code{\link{predict}}.

\code{\link[lme4]{sigma}}, from package \CRANpkg{lme4}.
}
\keyword{models}
