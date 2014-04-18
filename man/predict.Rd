\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{predict}
\alias{predict.hglm.fast}
\alias{sigma}
\alias{sigma.hglm.fast}
\title{
Prediction
}
\description{
\code{predict} gives empirical Bayes predictions of the response, while
\code{sigma} gives the dispersion parameter.
}
\usage{
\method{predict}{hglm.fast}(object, newdata = NULL, type = c("link", "response"),
          se.fit = FALSE, na.action = na.pass, ...)

\method{sigma}{hglm.fast}(object, ...)
}
\arguments{
  \item{object}{
    an \code{hglm.fast} object
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

The \code{sigma} function gives the dispersion parameter of the model.
}
\seealso{
\code{\link{predict}}.

\code{\link[lme4]{sigma}}, from package \CRANpkg{lme4}.
}
\keyword{models}
