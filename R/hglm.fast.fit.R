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


hglm.fast.fit <- function(x, y, group, weights = rep(1, nobs), start = NULL,
                          etastart = NULL, mustart = NULL,
                          offset = rep(0, nobs), family = gaussian(), 
                          control = list(), method = "firthglm.fit",
                          intercept = TRUE, standardize = TRUE, steps = 1)
{
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    nobs <- NROW(y)
    nvars <- ncol(x)

    # group-specific estimates
    m <- rdglm.group.fit(x = x, y = y, group = group, weights = weights,
                         start = start, etastart = etastart, mustart = mustart,
                         offset = offset, family = family, control = control,
                         method = method, intercept = intercept)
    ngroups <- m$ngroups

    # compute pooled dispersion estimates
    dispersion.tot <- dispersion.pooled(m$dispersion, m$df.residual)
    df.residual.tot <- sum(m$df.residual)
    dispersion <- rep(dispersion.tot, ngroups)

    if (dispersion.tot == 0)
        stop("cannot estimate dispersion (no residual degrees of freedom)")

    # compute group-sqecific precision square root (unpivoted)
    Rp <- as.list(rep(NULL, ngroups))
    for (i in seq_len(ngroups)) {
        qr.i <- m$qr[[i]]
        pivot.i <- qr.i$pivot
        rank.i <- qr.i$rank
        R.i <- qr.R(qr.i)
        Rp[[i]] <- matrix(0, rank.i, nvars)
        Rp[[i]][,pivot.i] <- R.i[seq_len(rank.i),]
    }

    # change coordinates so that average precision is identity
    if (standardize) {
        # compute averge precision of estimates
        prec.avg <- matrix(0, nvars, nvars)

        for (i in seq_len(ngroups)) {
            scale.i <- sqrt(1/dispersion[i])
            prec.avg <- prec.avg + crossprod(scale.i * Rp[[i]])
        }
        prec.avg <- prec.avg / ngroups

        # cov(coef) = (t(R) R)^{-1} = R^{-1} R^{-T}
        # cov(R x) = R cov(x) R^T
        # [cov(R x)]^{-1} = R^{-T} [cov(x)]^{-1} R^{-1}
        suppressWarnings({
            R <- chol(prec.avg, pivot = TRUE)
        })

        pivot <- attr(R, "pivot")
        rank <- attr(R, "rank")
        r1 <- seq_len(rank)
        R <- R[r1,r1,drop=FALSE]

    } else {
        R <- diag(nvars)
        pivot <- seq_len(nvars)
        rank <- nvars
        r1 <- seq_len(rank)
    }

    # compute standardized coeficients:
    #   coef[i,] <- R %*% m$coefficients[i,pivot[r1]]
    coef <- m$coefficients[,pivot[r1],drop=FALSE] %*% t(R)

    # compute group-specific subspace and precisions
    subspace <- as.list(rep(NULL, ngroups))
    precision <- as.list(rep(NULL, ngroups))
    for (i in seq_len(ngroups)) {
        prec.sqrt <- backsolve(R, t(Rp[[i]][,pivot[r1],drop=FALSE]),
                                transpose=TRUE)
        prec.sqrt.svd <- svd(prec.sqrt)

        subspace[[i]] <- prec.sqrt.svd$u
        precision[[i]] <- (prec.sqrt.svd$d)^2
    }

    # compute coefficient mean and covariance estimates
    est0 <- NULL
    suppressWarnings({
        for (s in seq_len(steps)) {
            est0 <- moment.est(coef, subspace, precision, dispersion, start.cov=est0$cov)
        }
    })
    est <- moment.est(coef, subspace, precision, dispersion, start.cov=est0$cov)
    mean <- est$mean
    mean.cov <- est$mean.cov
    cov <- est$cov

    # change back to original coordinates
    coef.mean <- rep(NA, nvars)
    coef.mean.cov <- matrix(NA, nvars, nvars)
    coef.cov <- matrix(NA, nvars, nvars)
    coef.mean[pivot[r1]] <- backsolve(R, mean)
    coef.mean.cov[pivot[r1],pivot[r1]] <- backsolve(R, t(backsolve(R, mean.cov)))
    coef.cov[pivot[r1],pivot[r1]] <- backsolve(R, t(backsolve(R, cov)))

    # set coordinate names
    names(coef.mean) <- xnames
    dimnames(coef.mean.cov) <- dimnames(coef.cov) <- list(xnames, xnames)

    z <- list(family = family, coefficient.mean = coef.mean,
              coefficient.mean.cov = coef.mean.cov, coefficient.cov = coef.cov,
              coefficients = coef, subspace = subspace, precision = precision,
              dispersion = dispersion.tot, df.residual = df.residual.tot,
              R = R, rank = rank, pivot = pivot, x = x, y = y, group = group,
              prior.weights = weights, offset = offset)
    z
}

