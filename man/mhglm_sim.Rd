\name{mhglm_sim}
\alias{mhglm_sim}
\title{
Simulate response patterns
}
\description{
Simulate response patterns for generalized linear models of \code{gaussian} or
\code{binomial} families, with both an intercept and slope covariate. Used
primarily for testing purposes. 
}
\usage{
mhglm_sim(n, m_per_level, sd_intercept, sd_slope,
          family = c("gaussian", "binomial"), seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
an integer scalar, the number of observations at the lowest grouping level.
}
  \item{m_per_level}{
an integer vector, the number of grouping levels nested under the level above.
}
  \item{sd_intercept}{
a numeric vector, the standard deviations of the intercept random effects. 
}
  \item{sd_slope}{
a numeric vector, the standard deviations of the slope random effects. 
}
  \item{family}{
a character scalar, either \code{"gaussian"} or \code{"binomial"}.
}
  \item{seed}{
a single value, interpreted as an integer, or \code{NULL} as in \code{set.seed}.
}
}
\details{
returns a data.frame with design matrix, response, and group levels. 
}
\examples{
mhglm_sim(n = 2, m_per_level = c(3, 3), sd_intercept = c(1, 2), 
        sd_slope = c(2, 1), family = "gaussian", seed = 123)

mhglm_sim(n = 2, m_per_level = c(3, 3), sd_intercept = c(1, 2), 
        sd_slope = c(2, 1), family = "binomial", seed = 123)
}
\keyword{datagen}
