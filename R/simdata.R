set.seed(12345)
n_level_1 <- 10
n_level_2 <- 5
n_level_3 <- 2
n_obs <- 30

x <- runif(n = n_level_1 * n_obs * n_level_2 * n_level_3, min = -1, max = 1)
g1 <- rep(seq_len(n_level_1), each = n_obs * n_level_2 * n_level_3)
g2 <- rep(seq_len(n_level_2 * n_level_1), each = n_obs * n_level_3)
g3 <- rep(seq_len(n_level_2 * n_level_1 * n_level_3), each = n_obs)

intercept_fix <- runif(1)
intercept_g1 <- rnorm(n_level_1)
intercept_g2 <- rnorm(n_level_1 * n_level_2, sd = 0.5)
intercept_g3 <- rnorm(n_level_1 * n_level_2 * n_level_3, sd = sqrt(0.25))

slope_fix <- runif(1)
slope_g1 <- rnorm(n_level_1)
slope_g2 <- rnorm(n_level_1 * n_level_2, sd = sqrt(0.5))
slope_g3 <- rnorm(n_level_1 * n_level_2 * n_level_3, sd = 0.25)

noise <- rnorm(length(x))

## The model is y ~ 1 + x + (1 + x | g1) + (1 + x | g2) + (1 + x |g3)
## nested groups, independent covariates.
y <- x * (slope_fix + slope_g1[g1] + slope_g2[g2] + slope_g3[g3]) +
  intercept_fix + intercept_g1[g1] + intercept_g2[g2] + intercept_g3[g3] + noise

simdata <- data.frame(x, g1, g2, g3, y)
# devtools::use_data(simdata)

#' Simulated data for testing \code{mhglm_ml}.
"simdata"
