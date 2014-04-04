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


moment.est <- function(coefficients, subspace, precision, dispersion,
                       start.cov = NULL)
{
    ngroups <- nrow(coefficients)
    dim <- ncol(coefficients)
    names <- colnames(coefficients)

    if (ngroups == 0L || dim == 0L){
        cov <- matrix(0, dim, dim)
        dimnames(cov) <- list(names, names)
        return(cov)
    }

    # compute mean estimate and covariance bias correction
    wt <- array(0, dim=c(ngroups, dim, dim))
    wtb <- matrix(0, ngroups, dim)
    bias <- array(0, dim=c(ngroups, dim, dim))

    for (i in seq_len(ngroups)) {
        u <- subspace[[i]]
        l <- precision[[i]]
        sigma2 <- dispersion[i]
        r <- length(l)

        if (r == 0L) {
            next
        }

        if (is.null(start.cov)) {
            w <- u %*% diag(l / (l + sigma2), r, r) %*% t(u)
            bias.i <- u %*% diag(sigma2 * l
                                 / (l + sigma2)^2, r, r) %*% t(u)
        } else {
            # implementation trick to avoid 1/li:
            # W = U (U^T Sigma U + a L^{-1})^{-1} U^T
            #   = U L^{1/2} (L^{1/2} U^T Sigma U L^{1/2} + a I)^{-1} L^{1/2} U^T
            #   = Us (Us^T Sigma Us + a I)^{-1} Us^T
            #
            # B = W U (a L^{-1}) U^T W
            #
            #   = U S (S U^T Sigma U S + a I)^{-1} S U^T
            #     U (a S^{-2}) U^T
            #     U S (S U^T Sigma U S + a I)^{-1} S U^T
            #
            #   =  a U S (S U^T Sigma U S + a I)^{-2} S U^T
            #
            s <- sqrt(l)
            us <- u %*% diag(s, r, r)
            cov.ii <- t(us) %*% start.cov %*% us
            h <- cov.ii + diag(sigma2, r, r)
            h.sqrt <- chol(h)
            hi.us.t <- backsolve(h.sqrt, t(us), transpose=TRUE)
            w <- t(hi.us.t) %*% hi.us.t

            h2i.us.t <- backsolve(h.sqrt, hi.us.t)
            bias.i <- sigma2 * t(h2i.us.t) %*% h2i.us.t
        }

        wt[i,,] <- w
        wtb[i,] <- w %*% coefficients[i,]
        bias[i,,] <- bias.i
    }

    wtot <- apply(wt, c(2,3), mean)
    mean <- pseudo.solve(wtot, colMeans(wtb))
    if (attr(mean, "deficient")) {
        warning("cannot solve mean moment equation due to rank deficiency")
    }
    mean.cov <- pseudo.solve(wtot) / ngroups
    attr(mean, "deficient") <- attr(mean.cov, "deficient") <- NULL

    # compute covariance estimate
    wt2 <- array(NA, dim=c(ngroups, dim^2, dim^2))
    wtb2 <- array(NA, dim=c(ngroups, dim, dim))

    for (i in seq_len(ngroups)) {
        w <- matrix(wt[i,,], dim, dim)
        wt2[i,,] <- kronecker(w, w)

        diff <- w %*% (coefficients[i,] - mean)
        wtb2[i,,] <- diff %*% t(diff)
    }

    wtot2 <- apply(wt2, c(2,3), mean)

    #wt.cov <- apply(wtb2 - bias, c(2, 3), mean)
    wt.cov <- apply(wtb2, c(2, 3), mean)
    wt.bias <- apply(bias, c(2, 3), mean)

    cov.vec <- pseudo.solve(wtot2, as.vector(wt.cov))
    bias.vec <- pseudo.solve(wtot2, as.vector(wt.bias))
    if (attr(cov.vec, "deficient")) {
        warning("cannot solve covariance moment equation due to rank deficiency")
    }
    cov <- matrix(cov.vec, dim, dim)
    bias <- matrix(bias.vec, dim, dim)

    cov <- 0.5 * (cov + t(cov)) # ensure symmetry
    bias <- 0.5 * (bias + t(bias))

    eigen.cov <- eigen(cov, symmetric=TRUE)
    l <- eigen.cov$values
    u <- eigen.cov$vectors[,l > 0, drop=FALSE]
    l <- l[l > 0]
    s <- sqrt(l)
    s.u.t <- t(u) * s
    sinv.u.t <- t(u) / s
    cov.bias <- sinv.u.t %*% bias %*% t(sinv.u.t)
    eigen.cov.bias <- eigen(cov.bias, symmetric=TRUE)
    l.bias <- eigen.cov.bias$values
    u.bias.t <- t(eigen.cov.bias$vectors) %*% s.u.t
    scale <- max(1, l.bias[1])
    cov.adj <- t(u.bias.t) %*% diag(1 - l.bias / scale, length(l.bias)) %*% u.bias.t

    cov <- proj.psd(cov.adj)  # ensure positive definite
    if (attr(cov, "modified") || length(l) < nrow(cov) || scale != 1)
        warning(paste("moment-based covariance matrix estimate is not positive"
                    , " semi-definite; using projection"
                    , sep=""))
    attr(cov, "modified") <- NULL

    names(mean) <- names
    dimnames(mean.cov) <- dimnames(cov) <- list(names, names)
    list(mean=mean, mean.cov=mean.cov, cov=cov)
}
