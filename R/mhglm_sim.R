mhglm_sim <- function(n, m_per_level, sd_intercept, sd_slope,
                    family = c("gaussian", "binomial"), seed) {
    set.seed(seed)
    total_n <- n * prod(m_per_level)
    n_levels <- length(m_per_level)

    groups <- Map(function(n, e) rep(seq_len(n), each = e),
                  cumprod(m_per_level),
                  n * rev(cumprod(rev(c(m_per_level[-1], 1)))))
    names(groups) <- paste0("g", seq_along(groups))

    fixed_effects <- setNames(rnorm(2), c("intercept", "slope"))
    ranef_intercept <- Map(rnorm, n = cumprod(m_per_level), sd = sd_intercept)
    ranef_slope <- Map(rnorm, n = cumprod(m_per_level), sd = sd_slope)

    total_ranef_intercept <- Reduce("+", Map("[", ranef_intercept, groups))
    total_ranef_slope <- Reduce("+", Map("[", ranef_slope, groups))

    x <- runif(n = total_n, min = -1, max = 1)
    z <- fixed_effects["intercept"] + total_ranef_intercept +
        x * (fixed_effects["slope"] + total_ranef_slope)

    y <- if (family == "gaussian") {
        z + rnorm(total_n)
    } else if (family == "binomial") {
        as.integer(runif(total_n) < plogis(z))
    } else {
        stop("family not supported")
    }

    cbind(x = x, as.data.frame(groups), y = y)
}
