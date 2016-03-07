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


# fit group-specific rank-deficient generalized linear models
rdglm.group.fit <- function(x, y, group, weights = rep(1, nobs), start = NULL,
                            etastart = NULL, mustart = NULL,
                            offset = rep(0, nobs), family = gaussian(),
                            control = list(), method = "firthglm.fit",
                            intercept = TRUE, parallel = FALSE)
{
    if(parallel) x <- bigmemory::as.big.matrix(x, shared=TRUE)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    nobs <- NROW(y)
    nvars <- ncol(x)

    group <- as.factor(group)
    levels <- levels(group)
    ngroups <- length(levels)

    coefficients <- matrix(NA, ngroups, nvars)
    dispersion <- rep(NA, ngroups)
    df.residual <- rep(NA, ngroups)

    rank <- rep(as.integer(NA), ngroups)
    qr <- as.list(rep(NULL, ngroups))

    group.int <- as.integer(group) # observation index => index of its group
    group.size <- tabulate(group.int, ngroups) # group index => # obs in that group
    subset <- lapply(group.size, integer) # group => vector of indices
    group.pos <- integer(ngroups)

    .Call(C_rdglm_index_loop, group.int, subset)

    if(parallel) {
        logging::loginfo("Fitting models in parallel", logger="mbest.mhglm.fit")
        results <- foreach(i=seq_len(ngroups)) %dopar% {
            j <- subset[[i]]
            yj <- if (is.matrix(y))
                y[j,,drop=FALSE]
            else y[j]

            model <- rdglm.fit(x = x[j,,drop=FALSE], y = yj,
                               weights = weights[j], start = start,
                               etastart = etastart[j], mustart = mustart[j],
                               offset = offset[j], family = family,
                               control = control, method = method,
                               intercept = intercept)

            return(list(model$coefficients,
                        model$dispersion,
                        model$df.residual,
                        model$rank,
                        model$qr))
        }

        coefficients <- t(vapply(results, function(x) x[[1]], numeric(nvars)))
        dispersion <- vapply(results, function(x) x[[2]], numeric(1))
        df.residual <- vapply(results, function(x) x[[3]], numeric(1))
        rank <- lapply(results, function(x) x[[4]])
        qr <- lapply(results, function(x) x[[5]])
        rm(results)
    } else {
        logging::loginfo("Fitting models in sequence", logger="mbest.mhglm.fit")
        for(i in seq_len(ngroups)) {
            j <- subset[[i]]
            yj <- if (is.matrix(y))
                y[j,,drop=FALSE]
            else y[j]

            model <- rdglm.fit(x = x[j,,drop=FALSE], y = yj,
                               weights = weights[j], start = start,
                               etastart = etastart[j], mustart = mustart[j],
                               offset = offset[j],
                               parallel=parallel, # cast x as matrix each time
                               family = family, control = control,
                               method = method, intercept = intercept)
            coefficients[i,] <- model$coefficients
            dispersion[i] <- model$dispersion
            df.residual[i] <- model$df.residual
            rank[[i]] <- model$rank
            qr[[i]] <- model$qr
        }
    }

    colnames(coefficients) <- xnames
    rownames(coefficients) <- levels
    names(dispersion) <- names(df.residual) <- names(rank) <- names(qr) <- levels

    list(ngroups = ngroups, group = group, coefficients = coefficients,
         dispersion = dispersion, df.residual = df.residual, rank = rank,
         qr = qr)
}

