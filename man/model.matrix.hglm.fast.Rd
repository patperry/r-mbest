\name{model.matrix.hglm.fast}
\alias{model.matrix.hglm.fast}
\alias{terms.hglm.fast}
\title{
Terms and Model Matrix
}
\description{
Get the terms or model matrix from an \code{hglm.fast} object.
}
\usage{
  \method{model.matrix}{hglm.fast}(object, type = c("fixed", "random"), ...)

  \method{terms}{hglm.fast}(x, type = c("fixed", "random"), ...)
}
\arguments{
  \item{object, x}{
    an \code{hglm.fast} object.
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
