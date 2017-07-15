\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{simdata}
\alias{simdata}
\title{
Simulate response patterns
}
\description{
Simulate response patterns for generalized linear models of \code{gaussian} or
\code{binomial} families, with both an intercept and slope covariate. Used
primarily for testing purposes. 
}
\usage{
simdata(n, m_per_level, sd_intercept, sd_slope, family = c("gaussian", "binomial"), seed)
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
simdata(n = 2, m_per_level = c(3, 3), sd_intercept = c(1, 2), 
        sd_slope = c(2, 1), family = "gaussian", seed = 123)
#             x g1 g2          y
# 1   0.5978497  1  1  3.1390310
# 2  -0.7562015  1  1  7.4841554
# 3   0.1218960  1  2  1.7800785
# 4  -0.5869372  1  2  1.9900380
# 5  -0.7449367  1  3 -0.5605606
# 6   0.5066157  1  3 -3.1459276
# 7   0.7900907  2  4 -0.2734820
# 8  -0.2510744  2  4 -3.7004708
# 9   0.3302304  2  5  1.8198573
# 10 -0.8103187  2  5 -2.7061211
# 11 -0.2320607  2  6  0.2972535
# 12 -0.4512327  2  6  0.5095146
# 13  0.6292801  3  7 -0.1551637
# 14 -0.1029673  3  7  1.0646329
# 15  0.6201287  3  8  0.3741042
# 16  0.6247790  3  8  0.7114451
# 17  0.5886846  3  9 -0.7806478
# 18 -0.1203366  3  9 -0.1418425

simdata(n = 2, m_per_level = c(3, 3), sd_intercept = c(1, 2), 
        sd_slope = c(2, 1), family = "binomial", seed = 123)
#             x g1 g2 y
# 1   0.5978497  1  1 1
# 2  -0.7562015  1  1 1
# 3   0.1218960  1  2 1
# 4  -0.5869372  1  2 1
# 5  -0.7449367  1  3 0
# 6   0.5066157  1  3 0
# 7   0.7900907  2  4 1
# 8  -0.2510744  2  4 0
# 9   0.3302304  2  5 1
# 10 -0.8103187  2  5 0
# 11 -0.2320607  2  6 1
# 12 -0.4512327  2  6 1
# 13  0.6292801  3  7 1
# 14 -0.1029673  3  7 0
# 15  0.6201287  3  8 1
# 16  0.6247790  3  8 1
# 17  0.5886846  3  9 0
# 18 -0.1203366  3  9 0

}
\keyword{datagen}
