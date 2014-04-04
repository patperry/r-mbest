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


ebayes.est <- function(x, offset = rep(0, nobs), family = gaussian(),
                       coefficient, subspace, precision, dispersion,
                       coefficient.mean, coefficient.cov)
{
    x <- as.matrix(x)
    nobs <- nrow(x)
    if (is.null(offset))
        offset <- rep(0, nobs)

    r <- length(precision)
    coef <- coefficient
    coef.mu <- coefficient.mean
    coef.cov <- coefficient.cov

    if (r == 0L) {
        coef.eb <- coefficient.mean
    } else {
        # implementation trick to avoid 1/li:
        # U (U^T Sigma U + a L^{-1})^{-1} U^T
        #   = U L^{1/2} (L^{1/2} U^T Sigma U L^{1/2} + a I)^{-1} L^{1/2} U^T
        #   = Us (Us^T Sigma Us + a I)^{-1} Us^T
        u <- subspace
        s <- sqrt(precision)
        us <- u %*% diag(s, r, r)
        cov.ii <- t(us) %*% coef.cov %*% us
        h <- cov.ii + diag(dispersion, r, r)
        w.diff <- us %*% solve(h, t(us) %*% (coef - coef.mu))
        coef.eb <- coef.mu + coef.cov %*% w.diff
    }

    eta <- offset + x %*% coef.eb
    mu <- family$linkinv(eta)

    list(coefficients = coef.eb, fitted.values = mu)
}


ebayes.group.est <- function(x, group, offset = rep(0, nobs), family = gaussian(),
                             coefficients, subspace, precision, dispersion,
                             coefficient.mean, coefficient.cov)
{
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    nobs <- nrow(x)
    nvars <- ncol(x)

    group <- as.factor(group)
    levels <- levels(group)
    ngroups <- length(levels)

    coefficients.eb <- matrix(NA, ngroups, nvars)
    colnames(coefficients.eb) <- xnames
    rownames(coefficients.eb) <- levels

    yhat <- numeric(nobs)

    group.int <- as.integer(group)
    group.size <- tabulate(group.int, ngroups)
    subset <- lapply(group.size, integer)
    group.pos <- integer(ngroups)
    for (o in seq_len(nobs)) {
        i <- group.int[o]
        j <- group.pos[i] + 1L
        subset[[i]][j] <- o
        group.pos[i] <- j
    }

    for (i in seq_len(ngroups)) {
        j <- subset[[i]]

        eb <- ebayes.est(x = x[j,,drop=FALSE], offset = offset[j], family = family,
                         coefficients[i,], subspace = subspace[[i]],
                         precision = precision[[i]], dispersion = dispersion[i],
                         coefficient.mean = coefficient.mean,
                         coefficient.cov = coefficient.cov)

        coefficients.eb[i,] <- eb$coefficients
        yhat[j] <- eb$fitted.values
    }

    list(coefficients = coefficients.eb, fitted.values = yhat)
}

