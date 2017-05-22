# Copyright 2014 Patrick O. Perry
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Compute 'sufficient statitics' for estimating fixed coef mean using moment method.
# Usually it will be used separately on each cores.
moment.est.mean.mapper <- function(coefficients, nfixed, subspace, precision,
                                   dispersion, start.cov = NULL)
{
    ngroups <- nrow(coefficients)
    dim <- ncol(coefficients)
    nrandom <- dim - nfixed

    fixed <- seq_len(nfixed)
    random  <- nfixed + seq_len(nrandom)

    if (is.null(start.cov))
        start.cov <- diag(nrandom)

    # compute mean estimate and covariance bias correction
    weight11.sum <- matrix(0,  nfixed, nfixed)
    weight1.coef.sum <- matrix(0,  nfixed, 1)

    for (i in seq_len(ngroups)) {
        u <- subspace[[i]]
        l <- precision[[i]]

        sigma2 <- dispersion[i]
        r <- length(l)

        if (r == 0L) {
            next
        }

        s <- sqrt(l)
        us <- u %*% diag(s, r, r)
        u1s <- us[fixed, , drop = FALSE]
        u2s <- us[random, , drop = FALSE]

        cov22 <- t(u2s) %*% start.cov %*% u2s
        w.inv <- cov22 + diag(sigma2, r, r)

        # pseudo.solve is more robust against 0 eigenvalues
        w11 <- u1s %*% pseudo.solve(w.inv, t(u1s))
        w12 <- u1s %*% pseudo.solve(w.inv, t(u2s))


        w1b <- (w11 %*% coefficients[i, fixed]
                + w12 %*% coefficients[i, random])

        stopifnot(max(abs(w1b)) <= 100)

        weight11.sum <- weight11.sum + w11
        weight1.coef.sum <- weight1.coef.sum + w1b
    }

    list(ngroups = ngroups, weight11.sum  = weight11.sum,
         weight1.coef.sum  = weight1.coef.sum )

}


# Estimate fixed effect's coefficient.
# Usually used as the 'combine' function in foreach function.
moment.est.mean.reducer <- function(ret, fixef.rank.warn)
{
    if (length(ret) == 0)
        stop("The input is empty.")

    ngroups <- 0
    wtot <- 0
    meanSum <- 0

    for (i in seq_len(length(ret))) {
        ngroups <- ngroups + ret[[i]]$ngroups
        wtot <- wtot + ret[[i]]$weight11.sum
        meanSum <- meanSum + ret[[i]]$weight1.coef.sum
    }
    wtot <- wtot / ngroups
    meanSum <- meanSum / ngroups

    mean <- pseudo.solve(wtot, meanSum)
    if (attr(mean, "deficient") & fixef.rank.warn) {
        warning(paste("cannot solve fixed effect moment equation due to rank",
                      "deficiency"))
    }
    mean.cov <- pseudo.solve(wtot) / ngroups
    attr(mean, "deficient") <- attr(mean.cov, "deficient") <- NULL

    list(mean = mean, mean.cov = mean.cov)
}

# Compute 'sufficient statitics' for estimating ranef cov using moment method.
# Usually it will be used separately on each cores.
moment.est.cov.mapper <- function(coefficients, nfixed, subspace, precision,
                                  dispersion, start.cov = NULL, diagcov = FALSE,
                                  mean)
{
    dim <- ncol(coefficients)
    ngroups <- nrow(coefficients)
    nrandom <- dim - nfixed
    fixed <- seq_len(nfixed)
    random  <- nfixed + seq_len(nrandom)

    if (is.null(start.cov))
        start.cov <- diag(nrandom)

    # compute mean estimate and covariance bias correction
    if (diagcov) {
        bias.sum <- matrix(0, nrandom, 1)
        weight22.2.sum <- matrix(0,  nrandom, nrandom)
        weight22.coef.2.sum <- matrix(0, nrandom, 1)
    } else {
        bias.sum <- matrix(0,  nrandom, nrandom)
        weight22.2.sum <- matrix(0,  nrandom^2, nrandom^2)
        weight22.coef.2.sum <- matrix(0,  nrandom, nrandom)
    }

    VWV12 <- list()
    VWV22 <- list()

    for (i in seq_len(ngroups)) {
        u <- subspace[[i]]
        l <- precision[[i]]

        sigma2 <- dispersion[i]
        r <- length(l)

        if (r == 0L) next

        s <- sqrt(l)
        us <- u %*% diag(s, r, r)
        u1s <- us[fixed, , drop = FALSE]
        u2s <- us[random, , drop = FALSE]

        cov22 <- t(u2s) %*% start.cov %*% u2s
        w.inv <- cov22 + diag(sigma2, r, r)

        # pseudo.solve is more robust against 0 eigenvalues
        w12 <- u1s %*% pseudo.solve(w.inv, t(u2s))
        w22 <- u2s %*% pseudo.solve(w.inv, t(u2s))

        R2.u2s.t <- pseudo.solve(w.inv, t(u2s))

        if (diagcov)
            B <- sigma2 * apply(R2.u2s.t ^ 2, 2, sum)
        else
            B <- sigma2 * t(R2.u2s.t) %*% R2.u2s.t

        bias.sum <- bias.sum + B

        if (diagcov) {
            diff <- (t(w12) %*% (coefficients[i, fixed] - mean)
                     + w22 %*% coefficients[i, random])

            weight22.2.sum <- weight22.2.sum + w22 ^ 2
            weight22.coef.2.sum  <- weight22.coef.2.sum + diff ^ 2
        } else {
            weight22.2.sum <- weight22.2.sum + kronecker(w22, w22)

            diff <- (t(w12) %*% (coefficients[i,fixed] - mean)
                     + w22 %*% coefficients[i,random])
            weight22.coef.2.sum <- weight22.coef.2.sum + diff %*% t(diff)
        }

        VWV12[[i]] <- w12
        VWV22[[i]] <- w22
    }

    list(ngroups = ngroups, weight22.2.sum = weight22.2.sum,
         weight22.coef.2.sum = weight22.coef.2.sum, bias.sum = bias.sum,
         VWV12 = VWV12, VWV22 = VWV22)
}


# Estimate ranef covariance matrix, if control$diagcov is TRUE, i.e. only estimate diagonal element in covariance matrix.
# Usually used as the 'combine' function in foreach function.
moment.est.cov.reducer <- function(ret, diagcov, cov.rank.warn, cov.psd.warn)
{
    if (length(ret) == 0)
        stop("The input is empty.")

    ngroups <- 0
    wtot2 <- 0
    wt.cov <- 0
    wt.bias <- 0

    for (i in seq_len(length(ret))) {
        ngroups <- ngroups + ret[[i]]$ngroups
        wtot2 <- wtot2 + ret[[i]]$weight22.2.sum
        wt.cov <- wt.cov + ret[[i]]$weight22.coef.2.sum
        wt.bias <- wt.bias + ret[[i]]$bias.sum
    }

    nrandom <- nrow(wt.bias)

    if (diagcov) {
        LHSinv <- pseudo.solve(wtot2)
        diag_1 <- LHSinv %*% wt.cov
        diag_2 <- LHSinv %*% wt.bias

        idx.use <- diag_2 !=0
        if ( any(idx.use)) {
            gamma <- min(min(diag_1[idx.use] / diag_2[idx.use]), 1)
        } else {
            gamma <- 0
        }

        cov.adj <- diag(c(diag_1 - gamma * diag_2), nrow = nrandom)
    } else {
        # construct an orthonormal basis for the space of symmetric
        # matrices
        q <- nrandom
        F <- matrix(0, q ^ 2, q * (q + 1) / 2)
        j <- 0
        for (k in seq_len(q)) {
            for (l in seq_len(k)) {
                j <- j + 1
                f <- matrix(0, q, q)
                if (k == l) {
                    f[k, l] <- 1
                } else {
                    f[k, l] <- 1 / sqrt(2)
                    f[l, k] <- 1 / sqrt(2)
                }
                F[, j] <- as.vector(f)
            }
        }

        # solve the moment equations
        tF.wtot2.F <- t(F) %*% wtot2 %*% F
        cov.vec <- pseudo.solve(tF.wtot2.F, t(F) %*% as.vector(wt.cov))
        bias.vec <- pseudo.solve(tF.wtot2.F, t(F) %*% as.vector(wt.bias))

        if (attr(cov.vec, "deficient") & cov.rank.warn) {
            warning(paste("cannot solve covariance moment equation due to rank",
                          "deficiency"))
        }

        # change back to original space
        cov <- matrix(F %*% cov.vec, nrandom, nrandom)
        bias <- matrix(F %*% bias.vec, nrandom, nrandom)

        # remove asymmetry arising from numerical errors
        cov <- 0.5 * (cov + t(cov))
        bias <- 0.5 * (bias + t(bias))

        eigen.cov <- eigen(cov, symmetric = TRUE)
        l <- eigen.cov$values
        u <- eigen.cov$vectors[, l > 0, drop = FALSE]
        l <- l[l > 0]

        if (length(l) == 0) {
            cov.adj <- diag(0, nrow = nrandom)
        } else {
            s <- sqrt(l)
            s.u.t <- t(u) * s
            sinv.u.t <- t(u) / s
            cov.bias <- sinv.u.t %*% bias %*% t(sinv.u.t)
            eigen.cov.bias <- eigen(cov.bias, symmetric = TRUE)
            l.bias <- eigen.cov.bias$values
            u.bias.t <- t(eigen.cov.bias$vectors) %*% s.u.t
            scale <- max(1, l.bias[1])
            cov.adj <- (t(u.bias.t)
                        %*% diag((scale - l.bias) / scale, length(l.bias))
                        %*% u.bias.t)
        }
    }

    cov <- proj.psd(cov.adj)  # ensure positive definite
    if (attr(cov, "modified") & cov.psd.warn)
        warning(paste("moment-based covariance matrix estimate is not positive",
                      "semi-definite; using projection"))
    attr(cov, "modified") <- NULL

    cov

}


moment.est <- function(coefficients, nfixed, subspace, precision, dispersion,
                       start.cov = NULL, parallel = FALSE, diagcov = FALSE,
                       fixef.rank.warn = FALSE, cov.rank.warn,
                       cov.psd.warn = TRUE)
{
    logging::loginfo("Estimating moments", logger = "mbest.mhglm.fit")
    ngroups <- nrow(coefficients)
    dim <- ncol(coefficients)
    names <- colnames(coefficients)
    if (ngroups == 0L || dim == 0L) {
        cov <- matrix(0, dim, dim)
        dimnames(cov) <- list(names, names)
        return(cov)
    }


    logging::loginfo("Computing mean estimate", logger = "mbest.mhglm.fit")
    i <- NULL

    # fixef
    if (parallel) {
        mean.info <- foreach(i = seq_len(ngroups)) %dopar% {
            moment.est.mean.mapper(coefficients[i, , drop = FALSE], nfixed,
                                   list(subspace[[i]]), list(precision[[i]]),
                                   dispersion[i], start.cov = start.cov)
        }

    } else {
        mean.info <- list(
            moment.est.mean.mapper(coefficients, nfixed, subspace, precision,
                                   dispersion, start.cov = start.cov))
    }

    est.mean <- moment.est.mean.reducer(mean.info, fixef.rank.warn)


    # ranef
    logging::loginfo("Computing covariance estimate",
                     logger = "mbest.mhglm.fit")
    if (parallel) {
        cov.info <- foreach(i = seq_len(ngroups)) %dopar% {
            moment.est.cov.mapper(coefficients[i, , drop = FALSE], nfixed,
                                  list(subspace[[i]]), list(precision[[i]]),
                                  dispersion[i],
                                  start.cov = start.cov, diagcov = diagcov,
                                  est.mean$mean)
            }
    } else {
        cov.info <- list(moment.est.cov.mapper(coefficients, nfixed,
                                               subspace, precision, dispersion,
                                               start.cov = start.cov,
                                               diagcov = diagcov,
                                               est.mean$mean))
    }

    est.cov <- moment.est.cov.reducer(cov.info, diagcov, cov.rank.warn,
                                      cov.psd.warn)

    list(mean = est.mean$mean, mean.cov = est.mean$mean.cov, cov = est.cov,
         mean.info = mean.info, cov.info = cov.info)
}
