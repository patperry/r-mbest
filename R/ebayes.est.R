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


ebayes.est <- function(coefficients, subspace, precision, dispersion,
                       coefficient.mean, coefficient.cov)
{
    coef <- coefficients
    coef.mu <- coefficient.mean
    coef.cov <- coefficient.cov
    r <- length(precision)

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

    coef.eb
}


ebayes.group.est <- function(coefficients, subspace, precision, dispersion,
                             coefficient.mean, coefficient.cov)
{
    ngroups <- nrow(coefficients)
    nvars <- ncol(coefficients)

    coefficients.eb <- matrix(NA, ngroups, nvars)
    dimnames(coefficients.eb) <- dimnames(coefficients)

    for (i in seq_len(ngroups)) {
        eb <- ebayes.est(coefficients[i,], subspace = subspace[[i]],
                         precision = precision[[i]], dispersion = dispersion[i],
                         coefficient.mean = coefficient.mean,
                         coefficient.cov = coefficient.cov)

        coefficients.eb[i,] <- eb
    }

    coefficients.eb
}

