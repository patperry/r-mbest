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


mhglm.fit <- function(x, z, y, group, weights = rep(1, nobs),
                      start = NULL, etastart = NULL, mustart = NULL,
                      offset = rep(0, nobs), family = gaussian(),
                      control = list(), intercept = TRUE, dispersion = NULL)
{
    control <- do.call("mhglm.control", control)
    x <- as.matrix(x)
    z <- as.matrix(z)
    xnames <- dimnames(x)[[2L]]
    znames <- dimnames(z)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    nobs <- NROW(y)
    nfixed <- ncol(x)
    nrandom <- ncol(z)
    nvars <- nfixed + nrandom
    fixed <- (seq_len(nvars) <= nfixed)
    random <- !fixed

    if (!is.null(start)) {
        if (length(start) != nfixed) {
            stop(gettextf(paste0("length of 'start' should equal %d",
                                 " and correspond to initial coefs for %s"),
                          nfixed, paste(deparse(xnames), collapse=", ")),
                 domain=NA)
        }
        start <- c(start, numeric(nrandom))
    }

    # group-specific estimates
    m <- rdglm.group.fit(x = cbind(x, z), y = y, group = group,
                         weights = weights, start = start,
                         etastart = etastart, mustart = mustart,
                         offset = offset, family = family,
                         parallel = control$parallel,
                         control = control$fit.control,
                         method = control$fit.method,
                         intercept = intercept)
    ngroups <- m$ngroups

    # compute pooled dispersion estimates
    if(is.null(dispersion)){
        dispersion.tot <- dispersion.pooled(m$dispersion, m$df.residual)
    } else {
        dispersion.tot <- dispersion
    }
    df.residual.tot <- sum(m$df.residual)
    dispersion <- rep(dispersion.tot, ngroups)

    #    if (dispersion.tot == 0)
    #        stop("cannot estimate dispersion (no residual degrees of freedom)")

    # compute group-sqecific precision square root (unpivoted)
    Rp <- as.list(rep(NULL, ngroups))

    if(control$parallel) {
        Rp <- foreach(i = seq_len(ngroups)) %dopar% {
            qr.i <- m$qr[[i]]
            return(unpivotRp(qr.i,nvars))
        }
    } else {
        for(i in seq_len(ngroups)) {
            qr.i <- m$qr[[i]]
            Rp[[i]] <- unpivotRp(qr.i,nvars)
        }
    }
    names(Rp) <- names(m$qr)

    # change coordinates so that average precision is identity
    if (control$standardize) {
        # compute averge precision of estimates

        if(control$diagcov){
            # if assuming independent covariates,
            # standardize each column separately.
            prec.avg <- rep(0, nvars)
            for (i in seq_len(ngroups)) {
                scale.i <- sqrt(1/dispersion[i])
                prec.avg <- prec.avg + apply((scale.i * Rp[[i]])^2, 2,sum)
            }
            prec.avg <- diag(prec.avg,nrow = nvars)

        } else {

            prec.avg <- matrix(0, nvars, nvars)
            if(control$parallel) {
                prec.avg <- foreach(i = seq_len(ngroups)) %dopar% {
                    scale.i <- sqrt(1/dispersion[i])
                    return(crossprod(scale.i * Rp[[i]]))
                }
                prec.avg <- Reduce('+', prec.avg)
            } else {
                for (i in seq_len(ngroups)) {
                    scale.i <- sqrt(1/dispersion[i])
                    prec.avg <- prec.avg + crossprod(scale.i * Rp[[i]])
                }
            }
        }

        prec.avg <- prec.avg / ngroups

        # cov(coef) = (t(R) R)^{-1} = R^{-1} R^{-T}
        # cov(R x) = R cov(x) R^T
        # [cov(R x)]^{-1} = R^{-T} [cov(x)]^{-1} R^{-1}
        suppressWarnings({
            R.fixed <- chol(prec.avg[fixed,fixed], pivot = TRUE)
            R.random <- chol(prec.avg[random,random], pivot = TRUE)
        })

        #pivot <- attr(R, "pivot")
        #rank <- attr(R, "rank")
        #r1.fixed <- seq_len(attr(R.fixed, "rank"))
        #r1.random <- seq_len(attr(R.random, "rank"))
        #R.fixed <- R.fixed[r1.fixed,r1.fixed,drop=FALSE]
        #R.random <- R.random[r1.random,r1.random,drop=FALSE]
    }
    else { # control.standardize==FALSE
        R.fixed <- diag(nfixed)
        attr(R.fixed, "pivot") <- seq_len(nfixed)
        attr(R.fixed, "rank") <- nfixed
        R.random <- diag(nrandom)
        attr(R.random, "pivot") <- seq_len(nrandom)
        attr(R.random, "rank") <- nrandom

        #pivot <- seq_len(nvars)
        #rank <- nvars
        #r1 <- seq_len(rank)
    }

    rank.fixed <- attr(R.fixed, "rank")
    rank.random <- attr(R.random, "rank")
    rank <- rank.fixed + rank.random

    r1.fixed <- seq_len(rank.fixed)
    r1.random <- seq_len(rank.random)
    r1 <- seq_len(rank)

    pivot.fixed <- which(fixed)[attr(R.fixed, "pivot")]
    pivot.random <- which(random)[attr(R.random, "pivot")]
    pivot <- c(pivot.fixed[r1.fixed],  pivot.random[r1.random],
               pivot.fixed[-r1.fixed], pivot.random[-r1.random])

    R <- matrix(0, rank, rank)
    R[r1.fixed, r1.fixed] <- R.fixed[r1.fixed, r1.fixed]
    (R[rank.fixed + r1.random, rank.fixed + r1.random]
     <- R.random[r1.random, r1.random])

    # compute standardized coeficients:
    #   coef[i,] <- R %*% m$coefficients[i,pivot[r1]]
    coef <- m$coefficients[,pivot[r1],drop=FALSE] %*% t(R)

    # compute group-specific subspace and precisions
    if(control$parallel) {
        results <- foreach(i = seq_len(ngroups)) %dopar% {
            if (nrow(Rp[[i]]) > 0L) {
                prec.sqrt <- backsolve(R, t(Rp[[i]][,pivot[r1],drop=FALSE]),
                                       transpose=TRUE)
                prec.sqrt.svd <- svd(prec.sqrt)
                return(list(subspace=prec.sqrt.svd$u,
                            precision=(prec.sqrt.svd$d)^2))
            } else {
                return(list(subspace=matrix(0, rank, 0), precision=numeric()))
            }
        }
        subspace <- lapply(results, function(x) x[['subspace']])
        precision <- lapply(results, function(x) x[['precision']])
        rm(results)
    } else {
        subspace <- as.list(rep(NULL, ngroups))
        precision <- as.list(rep(NULL, ngroups))
        for (i in seq_len(ngroups)) {
            if (nrow(Rp[[i]]) > 0L) {
                prec.sqrt <- backsolve(R, t(Rp[[i]][,pivot[r1],drop=FALSE]),
                                       transpose=TRUE)
                prec.sqrt.svd <- svd(prec.sqrt)
                subspace[[i]] <- prec.sqrt.svd$u
                precision[[i]] <- (prec.sqrt.svd$d)^2
            } else {
                subspace[[i]] <- matrix(0, rank, 0)
                precision[[i]] <- numeric()
            }
        }
    }
    names(subspace) <- names(Rp)
    names(precision) <- names(Rp)

    # compute coefficient mean and covariance estimates
    logging::loginfo("Computing mean and covariance estimates",
                     logger="mbest.mhglm.fit")
    est0 <- NULL
    suppressWarnings({
        for (s in seq_len(control$steps)) {
            est0 <- moment.est(coef, nfixed=rank.fixed, subspace, precision,
                               dispersion, start.cov=est0$cov,
                               parallel=control$parallel,
                               diagcov = control$diagcov,
                               fixef.rank.warn = control$fixef.rank.warn,
                               cov.rank.warn = control$cov.rank.warn,
                               cov.psd.warn = control$cov.psd.warn)
            logging::loginfo("Refining mean and covariance estimates",
                             logger="mbest.mhglm.fit")
        }
    })
    est <- moment.est(coef, nfixed=rank.fixed, subspace, precision, dispersion,
                      start.cov=est0$cov, parallel=control$parallel,
                      diagcov = control$diagcov,
                      fixef.rank.warn = control$fixef.rank.warn,
                      cov.rank.warn = control$cov.rank.warn,
                      cov.psd.warn = control$cov.psd.warn)

    mean <- est$mean
    mean.cov <- est$mean.cov
    cov <- est$cov

    # change back to original coordinates
    coef.mean <- rep(NA, nfixed)
    coef.mean.cov <- matrix(NA, nfixed, nfixed)
    coef.cov <- matrix(NA, nrandom, nrandom)
    (coef.mean[attr(R.fixed, "pivot")[r1.fixed]]
     <- backsolve(R.fixed, mean))
    (coef.mean.cov[attr(R.fixed, "pivot")[r1.fixed],
     attr(R.fixed, "pivot")[r1.fixed]]
    <- backsolve(R.fixed, t(backsolve(R.fixed, mean.cov))))
    (coef.cov[attr(R.random, "pivot")[r1.random],
     attr(R.random, "pivot")[r1.random]]
    <- backsolve(R.random, t(backsolve(R.random, cov))))

    # set coordinate names
    names(coef.mean) <- xnames
    dimnames(coef.mean.cov) <- list(xnames, xnames)
    dimnames(coef.cov) <- list(znames, znames)

    fit <- list(family = family, coefficient.mean = coef.mean,
                coefficient.mean.cov = coef.mean.cov,
                coefficient.cov = coef.cov, coefficients = coef,
                subspace = subspace, precision = precision,
                dispersion = dispersion.tot, df.residual = df.residual.tot,
                R = R, rank = rank, rank.fixed = rank.fixed,
                rank.random = rank.random, pivot = pivot, y = y, group = group,
                prior.weights = weights, offset = offset, nobs = nobs,
                mean.info = est$mean.info, cov.info = est$cov.info)

    fit
}

construct.reg <- function
## This is a unit function.
## This function takes in a list of mhglm fit objects,
## return a new data set (y,x,z) to be fit by the following
## hglm two-level model:
## y ~ x * beta + z * ranef_i + eps
## eps ~ N(0,I), i.i.d.
## ranef_i ~ N(0,Sigma), i.i.d
(fit.list, ##<< a list of objects returned by mhglm.fit
 nrandom ##<< number of random effects' covariates
 ){

    ngroups <- length(fit.list)

    # add a small positive number 1e-7 to singular values
    newdata <- lapply(fit.list,function(x){
                          omega <- x$mean.info[[1]]$weight11.sum
                          omegasvd <- svd(omega)
                          omegasvd$d <- pmax(1e-7, omegasvd$d)
                          omega_sqrt <- omegasvd$v %*% diag(sqrt(omegasvd$d), nrow = nrow(omega)) %*% t(omegasvd$v)

                          newy <- omega_sqrt %*% x$coefficient.mean
                          newxz <- omega_sqrt

                          colnames(newxz) <- names(x$coefficient.mean)
                          list(newy = newy, newxz = newxz) })

    nfixed <- length(newdata[[1]]$newy) - nrandom

    # concatenate all ys
    newy <- Reduce("c",lapply(newdata, function(x) x$newy))

    # first nfixed columns are 'fixed effects covariates'
    newx <- Reduce("rbind",lapply(newdata, function(x) x$newxz[,(1:nfixed),drop = FALSE]))

    # last nrandom columns are 'random effects covariates'
    newz <- Reduce("rbind",lapply(newdata, function(x) x$newxz[,(nfixed + (1:nrandom)),drop = FALSE]))

    # generate group indicator and turn into factors
    levels <- names(fit.list)
    newgroup <- gl(ngroups, (nfixed + nrandom),labels = levels)

    return(list(y = newy, x = newx, z = newz,group = newgroup))
}


fit.recursive <- function
## This is a recursive function. Going from bottom to top.
## This function takes in a tree-structured list of mhglm fit objects,
## - 1. go to the nodes where its children has mhglm fit objects,
## - 2. construct a new regression problem by collecting information from its children,
## - 3. fit mhglm to the newly constructed regression problem.
## One call to the function will only go one level up.
## In order to process n-levels hierachical structure, one need to call this n-1 times.
## Example: suppose the nested group is g1/g2/g3, then this is a 3-levels model and we
## call the function 2 times.
(fit.tree,##<< a tree-structured list of mhglm fit objects
 nrandom, ##<< number of random effects' covariates
 control ##<< control parameter from the main function call
 ){
    if(class(fit.tree) == 'mhglmfit'){
        return(fit.tree)
    } else if(class(fit.tree[[1]]) != 'mhglmfit'){
        # Need to go one level down
        fit.tree <- lapply(fit.tree, function(x) fit.recursive(x, nrandom, control))
    } else {
        # Reach the right level.
        # construct new regression problem
        newdata <- construct.reg(fit.tree,nrandom)

        # Fit standard mhglm to it.
        # Note
        # - family = gaussian(), regardless of what family the main model is.
        # - set dispersion to a small constant to ensure positive semi-definite.
        ret <- mhglm.fit(x = newdata$x, z = newdata$z,
                         y = newdata$y, group = newdata$group,
                         family = gaussian(), control = control,
                         dispersion = 1)
        fit.tree <- c(fit.tree,ret)

        # Set the flag so next time the function is called,
        # it will work on the parent node.
        class(fit.tree) <- 'mhglmfit'
    }
    return(fit.tree)
}

avg.coef.cov <- function
## This is a recursive function.
## This function takes in a tree-structured list of mhglm fit objects,
## returns the estimated ranef covariance for the given level.
## It weights the random effects covariance estimate by number of subgroups,
## and sum over all the nodes on that level.
(fit.tree, ##<< a tree-structured list of mhglm fit objects
 r ##<< for which level to compute the ranef covariance estimate.
 ### 1 <= r <= n-1, where n is the number of levels
 ){

    if(r==1){
        ngroups <- length(fit.tree$subspace)
        coefficient.cov.weighted <- ngroups * fit.tree$coefficient.cov
        ret <- list(ngroups = ngroups, coefficient.cov.weighted = coefficient.cov.weighted)
        return(ret)
    } else if (r>1) {
        ret<- lapply(Filter(function(x) class(x) == 'mhglmfit', fit.tree),
                     function(x) avg.coef.cov(x,r-1))
        ngroups <- Reduce("+", lapply(ret,function(x) x$ngroups))
        coefficient.cov.weighted <- Reduce("+", lapply(ret,function(x) x$coefficient.cov.weighted))
        ret2 <- list(ngroups = ngroups, coefficient.cov.weighted = coefficient.cov.weighted)
        return(ret2)
    } else {
        stop("Unvalid r.")
    }
}


avg.dispersion <- function
## This function is similar to avg.coef.cov, but to compute dispersion.
## Also it will only be used to computes dispersion for the bottom level.
(fit.tree, ##<< a tree-structured list of mhglm fit objects
 r ##<< for which level to compute the ranef covariance estimate.
 ### 1 <= r <= n-1, where n is the number of levels
 ){
    if(r==1){
        return(c(fit.tree$df.residual, fit.tree$df.residual * fit.tree$dispersion))
    } else if (r>1) {
        ret <- lapply(Filter(function(x) class(x) == 'mhglmfit', fit.tree),
                      function(x) avg.dispersion(x,r-1))

        df.residual.tot <- Reduce("+", lapply(ret,function(x) x[1]))
        dispersion.tot <- Reduce("+", lapply(ret,function(x) x[2]))
        return(c(df.residual.tot, dispersion.tot))
    } else {
        stop("Unvalid r.")
    }
}


mhglm.fit.bottom <- function
## This is a recursive function.
## It takes in data and fit mhglm to the lowest level of grouping.
## It returns a tree-structured list of mhglm fit objects.
(x, z, y, group, weights = rep(1, nobs),
 start = NULL, etastart = NULL, mustart = NULL,
 offset = rep(0, nobs), family = gaussian(),
 control = list(), intercept = TRUE
 ){

    control <- do.call("mhglm.control", control)
    x <- as.matrix(x)
    z <- lapply(z, as.matrix)
    ngroupslevel <- length(group)

    if(ngroupslevel ==1){
        m <- mhglm.fit(x = x, z = z[[1]], y = y, group = group[[1]], weights = weights,
                       start = start, etastart = etastart, mustart = mustart,
                       offset = offset, family = family,
                       control = control, intercept = intercept)
        class(m) <- 'mhglmfit'
        return(m)

    } else {

        ngroups <- nlevels(group[[1]])
        levels <- levels(group[[1]])
        fit.bottom <- as.list(rep(NULL, ngroups))
        subsets <- .Call(C_group_subsets, group[[1]], ngroups) # group => indices

        for(i in seq_len(ngroups)) {
            j <- subsets[[i]]
            yj <- if (is.matrix(y)) y[j,,drop=FALSE] else y[j]
            xj <- cbind(x[j,,drop = FALSE],
                        (z[[1]])[j,,drop = FALSE])

            groupj <- list()
            zj <- list()
            for(r in seq_len(ngroupslevel-1)){
                groupj[[r]] <- droplevels(group[[r+1]][j])
                zj[[r]] <- z[[r+1]][j,,drop = FALSE]
            }

            m <- mhglm.fit.bottom(x=xj, z = zj, y = yj,
                                  group = groupj, weights = weights[j],
                                  start = start, etastart = etastart, mustart = mustart,
                                  offset = offset[j], family = family,
                                  control = control, intercept = intercept)
            fit.bottom[[i]]<- m
        }

        names(fit.bottom) <- levels
        return(fit.bottom)

    }
}

mhglm_ml.fit <- function
#mhglm.fit.multilevel <- function
## This is the main function.
## This function fit multilevel GLM with more than 2 levels.
## TODO: can it fit two-level model?
## All the arguments except 'group' is same as original function.
## The 'group' argument must be a data frame (not necessarily consists
## of factors). But it will be turned into data frame of factors for later use.
(x, z, y, group, weights = rep(1, nobs),
 start = NULL, etastart = NULL, mustart = NULL,
 offset = rep(0, nobs), family = gaussian(),
 control = list(), intercept = TRUE
 ){

    control <- do.call("mhglm.control", control)
    x <- as.matrix(x)
    z <- lapply(z, as.matrix)
    xnames <- dimnames(x)[[2L]]
    znames <- lapply(z,function(x){dimnames(x)[[2L]]})

    nobs <- NROW(y)
    nfixed <- ncol(x)
    nrandom <- lapply(z,ncol)
    ngroupslevel <- if(is.list(z)){length(z)} else 0

    #----------
    # Right now we do not support control$standardize = TRUE
    if (control$standardize == TRUE){
        print('control$standardize = TRUE is currently not allowed. Set standardize to FALSE')
        control$standardize = FALSE
    }


    #----------
    # 'group' input needs special care.
    # Turn it into a list of factors
    if(!is.list(group)){
        stop('group must be a data.frame')
    } else{
        group <- lapply(group,factor)
    }


    #----------
    # run black box algorithm on groups
    fit.bottom <- mhglm.fit.bottom(x = x, z = z, y = y, group = group, weights = weights,
                                   start = start, etastart = etastart, mustart = mustart,
                                   offset = offset, family = family,
                                   control = control, intercept = intercept)

    #----------
    # recursively regress group specific estimate on true fixef and raneef
    r <- ngroupslevel
    while ( r>1){
        fit.bottom <- fit.recursive(fit.bottom, nrandom[[r-1]],control)
        r <- r-1
    }

    #----------
    # pool dispersion together
    dispersion.info <- avg.dispersion(fit.bottom,ngroupslevel)
    dispersion <- dispersion.info[2]/dispersion.info[1]

    #----------
    # Take average of estimated coefficient.cov
    coef.cov <- as.list(rep(NULL,ngroupslevel))
    for(r in seq_len(ngroupslevel)){
        coef.cov.info <- avg.coef.cov(fit.bottom,r)
        coef.cov[[r]] <- coef.cov.info[[2]]/coef.cov.info[[1]]
        dimnames(coef.cov[[r]]) <- list(znames[[r]], znames[[r]])
    }


    #----------
    # Set coordinate names
    names(fit.bottom$coefficient.mean) <- xnames
    dimnames(fit.bottom$coefficient.mean.cov) <- list(xnames, xnames)

    fit <- list(family = family,
                df.residual = dispersion.info[1],
                dispersion = dispersion,
                prior.weights = weights,
                fit = fit.bottom,
                coef.cov.all = coef.cov)

    fit
}




