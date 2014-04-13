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


ebayes.est <- function(coefficients, nfixed, subspace, precision, dispersion,
                       coefficient.mean, coefficient.cov)
{
    coef <- coefficients
    coef.mu <- coefficient.mean
    coef.cov <- coefficient.cov
    r <- length(precision)
    nrandom <- nrow(coef.cov)

    fixed <- seq_len(nfixed)
    random <- nfixed + seq_len(nrandom)

    if (r == 0L) {
        coef.eb <- numeric(nrandom)
    } else {
        # implementation trick to avoid 1/li:
        # U (U^T Sigma U + a L^{-1})^{-1} U^T
        #   = U L^{1/2} (L^{1/2} U^T Sigma U L^{1/2} + a I)^{-1} L^{1/2} U^T
        #   = Us (Us^T Sigma Us + a I)^{-1} Us^T
        u <- subspace
        s <- sqrt(precision)
        us <- u %*% diag(s, r, r)
        u1s <- us[fixed,,drop=FALSE]
        u2s <- us[random,,drop=FALSE]

        cov.ii <- t(u2s) %*% coef.cov %*% u2s
        h <- cov.ii + diag(dispersion, r, r)
        w.diff <- u2s %*% solve(h, (t(u1s) %*% (coef[fixed] - coef.mu)
                                    + t(u2s) %*% coef[random]))
        coef.eb <- coef.cov %*% w.diff
    }

    coef.eb
}


ebayes.group.est <- function(coefficients, nfixed, subspace, precision, dispersion,
                             coefficient.mean, coefficient.cov)
{
    ngroups <- nrow(coefficients)
    nvars <- ncol(coefficients) - nfixed

    coefficients.eb <- matrix(NA, ngroups, nvars)
    dimnames(coefficients.eb) <- dimnames(coefficients)

    for (i in seq_len(ngroups)) {
        eb <- ebayes.est(coefficients[i,],
                         nfixed = nfixed,
                         subspace = subspace[[i]],
                         precision = precision[[i]], dispersion = dispersion[i],
                         coefficient.mean = coefficient.mean,
                         coefficient.cov = coefficient.cov)

        coefficients.eb[i,] <- eb
    }

    coefficients.eb
}

