\name{model.matrix.mhglm}
\alias{model.matrix.mhglm}
\alias{terms.mhglm}
\title{
Terms and Model Matrix
}
\description{
Get the terms or model matrix from an \code{mhglm} object.
}
\usage{
  \method{model.matrix}{mhglm}(object, type = c("fixed", "random"), ...)

  \method{terms}{mhglm}(x, type = c("fixed", "random"), ...)
}
\arguments{
  \item{object, x}{
    an \code{mhglm} object.
  }
  \item{type}{
    which terms to get (for the fixed or for the random effects).
  }
  \item{\dots}{
    further arguments passed to or from other methods.
  }
}
\seealso{
\code{\link{model.matrix}}, \code{\link{terms}}
}
\keyword{models}
