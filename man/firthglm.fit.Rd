\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{firthglm.fit}
\alias{firthglm.fit}
\alias{firthglm.control}
\title{
Fitting Generalized Linear Models with Firth's Bias Reduction
}
\description{
A drop-in replacement for \code{\link{glm.fit}} which uses Firth's
bias-reduced estimates instead of maximum likelihood.
}
\usage{
firthglm.fit(x, y, weights = rep(1, nobs),
             start = NULL, etastart = NULL, mustart = NULL,
             offset = rep(0, nobs), family = gaussian(),
             control = list(...), intercept = TRUE, singular.ok = TRUE, ...)

firthglm.control(epsilon = 1e-8, maxit = 25, qr.tol = 1e-7,
                 improve.tol = 1e-4, curvature.tol = 0.9,
                 linesearch.method = "linesearch",
                 linesearch.maxit = 20, trace = FALSE)
}
\arguments{
  \item{x, y, weights, start, etastart, mustart, offset, family, control, intercept, singular.ok, \dots}{
    arguments that have the same functions as for \code{\link{glm.fit}}.
  }
  \item{qr.tol}{
    tolerance parameter for determining the rank of \code{x}.
  }
  \item{epsilon, maxit}{
    convergence parameters for the quasi-Newton method.
  }
  \item{linesearch.method}{
    line search methods (one of "linesearch", "backtrack", or "blindsearch")
  }
  \item{improve.tol, curvature.tol, linesearch.maxit}{
    tolerance parameters for the linesearch procedure.
  }
  \item{trace}{
    logical indicating if output should be produced for each iteration.
  }
}
\details{
Firth's modified score function gives rise to estimates with smaller biases
than their maximum likelihood counterparts.  Unlike the maximum likelihood
estimates, if the design matrix is of full rank, then the Firth bias-reduced
estimate is finite.

By default, the fitting procedure uses a quasi-Newton optimization method,
with a More-Thuente linesearch.
}
\value{
\code{firthglm.fit} returns an object having the same components that a call
to \code{\link{glm.fit}} would produce.
}
\references{
  Firth, D. (1993) Bias reduction of maximum likelihood estimates.
  \emph{Biometrika} \bold{80}, 27-–38.

  More, J. J. and Thuente, D. J. (1994) Line search algorithms with
  guaranteed sufficient decrease.
  \emph{ACM Transactions on Mathematical Software} \bold{20} 286--307. 
}
\author{
Patrick O. Perry
}
\note{
Currently, only families with canonical link functions are supported.
}
\seealso{
\code{\link[logistf]{logistf}} (package \CRANpkg{logistf}) and
\code{\link[brglm]{brglm}} (package \CRANpkg{brglm}) for alternative
implementations of Firth's bias-reduced estimators.
}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

## Use firthglm.fit instead of glm.fit:
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),
               method="firthglm.fit")
summary(glm.D93)
}
\keyword{models}
\keyword{regression}
