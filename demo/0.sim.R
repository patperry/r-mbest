library(devtools)
library(MASS)

gen.replicate <- function(ngroup, nobs, nobs.dispersion, nfixef, fixef.mean.df,
                          fixef.scale, nranef, ranef.cov.df, ranef.scale,
                          dispersion, family) {
    # we need to be strict about integer vs. numeric for the md5sum to match
    ngroup <- as.integer(ngroup)
    nobs <- as.integer(nobs)
    nobs.dispersion <- as.numeric(nobs.dispersion)
    nfixef <- as.integer(nfixef)
    fixef.mean.df <- as.numeric(fixef.mean.df)
    fixef.scale <- as.numeric(fixef.scale)
    nranef <- as.integer(nranef)
    ranef.cov.df <- as.numeric(ranef.cov.df)
    ranef.scale <- as.numeric(ranef.scale)
    dispersion <- as.numeric(dispersion)

    # family
    if (is.character(family)) {
        family <- get(family, mode = "function", envir = parent.frame())()
    } else if (is.function(family)) {
        family <- family()
    } else if (is.null(family$family)) {
        stop("'family' not recognized")
    }

    if (family$family %in% c("binomial", "poisson") && dispersion != 1) {
        warning(sprintf("dispersion parameter is ignored for '%s' family",
                        family$family))
    }

    # generate coefficient mean from t
    if (is.na(fixef.mean.df)) {
        fixef <- numeric(nfixef)
    } else {
        fixef <- fixef.scale * rt(nfixef, df = fixef.mean.df)
    }

    # generate global coefficient covariance from inverse Wishart
    if (is.na(ranef.cov.df)) {
        ranef.cov <- diag(ranef.scale ^ 2, nranef)
        ranef.cov.sqrt <- diag(ranef.scale, nranef)
    } else {
        rand_wish <- solve(rWishart(1, ranef.cov.df, diag(nranef))[, , 1])
        ranef.cov <- diag(diag(ranef.scale ^ 2 * rand_wish), nrow = nranef)
        ranef.cov.sqrt <- chol(ranef.cov)
    }


    # generate coefficients
    u <- matrix(rnorm(ngroup * nranef), ngroup, nranef)
    ranef <- u %*% ranef.cov.sqrt

    # generate sampling rates
    rate.mean <- nobs / ngroup
    if (nobs.dispersion == Inf) {
        sample.rate <- rep(rate.mean, ngroup)
    } else {
        sample.rate <- rgamma(ngroup, nobs.dispersion,
                              nobs.dispersion / rate.mean)
        stopifnot(all(sample.rate > 0))
    }

    # generate group
    suppressWarnings({ # ignore warning about using Walker's alias method
      group <- sample.int(ngroup, nobs, replace = TRUE, prob = sample.rate)
    })


    # generate feature  matrices with Pr(x[i,j] = +1) = P(x[i,j] = -1) = 1/2,
    z <- matrix(sample(c(-1, +1), nobs * nranef, replace = TRUE), nobs, nranef)
    x <- z

    # compute linear predictors and generate observations
    eta <- drop(x %*% fixef) + rowSums(z * ranef[group, ])
    mu <- family$linkinv(eta)

    if (family$family == "gaussian") {
        y <- rnorm(nobs, mean = mu, sd = sqrt(dispersion))
    } else if (family$family == "binomial") {
        y <- as.numeric(rbinom(nobs, 1, mu))
    } else if (family$family == "poisson") {
        y <- rpois(nobs, mu)
    } else {
        stop(sprintf("family '%s' not supported", family$family))
    }

    list(ngroup = ngroup, nobs = nobs,
         fixef = fixef, ranef = ranef, ranef.cov.sqrt = ranef.cov.sqrt,
         dispersion = dispersion, sample.rate = sample.rate,
         group = group, x = x, z = z, y.mean = mu, y = y, family = family)
}

set.seed(0)
r <- gen.replicate(ngroup = 10,  nobs = 10000, nobs.dispersion = 1, nfixef = 5,
                   fixef.mean.df = 4, fixef.scale = 1, nranef = 5,
                   ranef.cov.df = 10, ranef.scale = 1, dispersion = 1,
                   family = family)


group <- as.factor(r$group)
ngroups <- nlevels(group)
levels <- levels(group)
subsets <- .Call(C_group_subsets, group, ngroups) # group => indices

for (i in seq_len(length(subsets))) {
  idx <- subsets[[i]]
  data_filename <- paste0(datafile_folder, "data.", i, ".rds", sep = "")
  saveRDS(list(x = r$x[idx, ], z = r$z[idx, ], y = r$y[idx],
               group = factor(r$group[idx])), file = data_filename)
}

nfiles <- length(subsets)
