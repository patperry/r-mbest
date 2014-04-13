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


hglm.fast <- function(formula, group, family = gaussian, data,
                      weights, subset, na.action, start = NULL, etastart,
                      mustart, offset, control = list(), contrasts = NULL,
                      method = "firthglm.fit", standardize = TRUE, steps = 1,
                      x = FALSE, y = TRUE)
{
    # call
    call <- match.call()

    # family
    if (is.character(family))
        family <- get(family, mode="function", envir=parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }

    # data
    if (missing(data))
        data <- environment(formula)

    # model frame
    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula", "group", "data", "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    # method
    if (!is.character(method) && !is.function(method))
        stop("invalid 'method' argument")
    if (identical(method, "firthglm.fit"))
        control <- do.call("firthglm.control", control)
    if (identical(method, "glm.fit"))
        control <- do.call("glm.control", control)

    # response
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }

    # design matrix
    mt <- attr(mf, "terms")
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)

    # group
    group <- as.factor(model.extract(mf, "group"))
    if (length(group) != NROW(Y))
        stop(gettextf("number of groups is %d should equal %d (number of observations)", 
            length(group), NROW(Y)), domain = NA)

    # weights
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")

    # offset
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }

    # starting values
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")

    # group-specific estimates
    z <- hglm.fast.fit(x.fixed = X, x.random = X, y = Y, group = group, weights = weights,
                       start = start, etastart = etastart, mustart = mustart,
                       offset = offset, family = family, control = control,
                       method = method, intercept = attr(mt, "intercept") > 0L,
                       standardize = standardize, steps = steps)
    z$call <- call
    if (x)
        z$x <- X
    if (!y)
        z$y <- NULL
    class(z) <- "hglm"

    z
}


fixef.hglm <- function(object, ...)
{
    object$coefficient.mean
}

vcov.hglm <- function(object, ...)
{
    object$coefficient.mean.cov
}

VarCorr.hglm <- function(x, sigma = 1, rdig = 3)
{
    vc <- x$coefficient.cov
    stddev <- sqrt(diag(vc))
    cor <- scale(vc, center=FALSE, scale=stddev) / stddev

    attr(cor, "scaled:scale") <- NULL
    attr(vc, "stddev") <- stddev
    attr(vc, "correlation") <- cor

    varcor <- list()
    varcor[[as.character(x$call[["group"]])]] <- vc
    attr(varcor, "sc") <- sigma * sqrt(x$dispersion)
    attr(varcor, "useSc") <- !(x$family$family %in% c("binomial", "poisson"))
    class(varcor) <- "VarCorr.hglm"
    varcor
}

ranef.hglm <- function(object, condVar = FALSE, ...)
{
    if (condVar)
        warning("condVar is not implemented")

    nvars <- ncol(object$coefficients)
    xnames <- names(object$coefficient.mean)
    ngroups <- nrow(object$coefficients)
    gnames <- rownames(object$coefficients)

    R <- object$R
    pivot <- object$pivot
    rank <- object$rank
    r1 <- seq_len(rank)
    nfixed <- length(object$coefficient.mean)
    nrandom <- nvars - nfixed
    fixed <- seq_len(nfixed)
    random <- nfixed + seq_len(nrandom)

    R.fixed <- R[fixed,fixed,drop=FALSE]
    pivot.fixed <- pivot[fixed]
    coef.mean1 <- drop(R.fixed %*% object$coefficient.mean[pivot.fixed])

    R.random <- R[random,random,drop=FALSE]
    pivot.random <- pivot[random] - nfixed
    coef.cov1 <- (R.random
                  %*% object$coefficient.cov[pivot.random,
                                             pivot.random,
                                             drop=FALSE] %*% t(R.random))

    coef1 <- ebayes.group.est(coefficients=object$coefficients,
                              nfixed=nfixed,
                              subspace=object$subspace,
                              precision=object$precision,
                              dispersion=rep(object$dispersion, ngroups),
                              coefficient.mean=coef.mean1,
                              coefficient.cov=coef.cov1)

    # change back to original coordinates
    coef <- matrix(NA, ngroups, nrandom)
    coef[,pivot.random] <- t(backsolve(R.random, t(coef1)))
    colnames(coef) <- colnames(object$coefficient.cov)
    rownames(coef) <- gnames

    re <- list()
    re[[as.character(object$call[["group"]])]] <- coef
    class(re) <- "ranef.hglm"
    re
}

summary.hglm <- function(object, ...)
{
    # fixed effects
    vcov <- vcov(object)
    coefs <- cbind("Estimate" = fixef(object),
                   "Std. Error" = sqrt(diag(vcov)))
    coefs <- cbind(coefs, (cf3 <- coefs[,1]/coefs[,2]), deparse.level=0)
    colnames(coefs)[3] <- "z value"
    coefs <- cbind(coefs, "Pr(>|z|)" =
                          2*pnorm(abs(coefs[,3]), lower.tail = FALSE))

    # random effects
    varcor <- VarCorr(object)

    structure(list(call = object$call, family = object$family,
                   coefficients = coefs, dispersion = model$dispersion,
                   vcov = vcov, varcor = varcor),
              class = "summary.hglm")
}

print.VarCorr.hglm <- function (x, digits = max(3, getOption("digits") - 2),
                                var.print = FALSE, ...)
{
    dims <- sapply(x, ncol)
    pmax <- max(1L, max(dims))

    table <- matrix("", sum(dims + attr(x, "useSc")), pmax + 3L)

    rownames(table) <- rep("", nrow(table))
    colnames(table) <- c("Groups", "Name", "Variance", "Std.Dev", rep("", pmax - 1L))
    if (pmax > 1L) {
        colnames(table)[5L] <- "Corr"
    }

    off <- 0L
    for (i in seq_along(x)) {
        xx <- x[[i]]
        stddev <- attr(xx, "stddev")
        cor <- attr(xx, "correlation")
        p <- ncol(cor)

        tab <- matrix("", p, p + 3L)
        tab[1L,1L] <- names(x)[[i]]
        tab[,2L] <- names(stddev)
        tab[,3L] <- format(diag(xx), digits = digits, ...)
        tab[,4L] <- format(stddev, digits = digits, ...)

        if (p > 1L) {
            cor.str <- format(cor, digits = max(2L, digits - 2L), ...)
            cor.str[row(cor.str) <= col(cor.str)] <- ""
            tab[,5L:(p+3L)] <- cor.str[,-p]
        }
        table[off + seq_len(p), seq_len(p + 3L)] <- tab
        off <- off + p
    }

    if (attr(x, "useSc")) {
        sigma <- attr(x, "sc")
        table[off + 1L, 1L] <- "Residual"
        table[off + 1L, 3L] <- format(sigma^2, digits = digits, ...)
        table[off + 1L, 4L] <- format(sigma, digits = digits, ...)
    }

    if (!var.print) {
        table <- table[,-3L]
    }

    print.table(table)
}



print.summary.hglm <- function(x, digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

    cat("\nRandom effects:\n")
    print(x$varcor, var.print=TRUE)

    cat("\nFixed effects:\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)

    cat("\n")
}





print.hglm <- function(x, digits = max(3L, getOption("digits") - 3L),
                       signif.stars = getOption("show.signif.stars"), ...)
{
    names <- names(x$coefficient.mean)

    # random effects
    var <- pmax(0, diag(x$coefficient.cov))
    sd <- sqrt(var)
    random <- cbind(var, sd)
    colnames(random) <- c("Variance", "Std. Dev.")
    rownames(random) <- names

    # fixed effects
    est <- x$coefficient.mean
    est.se <- sqrt(pmax(0, diag(x$coefficient.mean.cov)))
    tstat <- est / est.se
    if (x$family$family %in% c("binomial", "poisson")) {
        pval <- 2 * pnorm(-abs(tstat))
    } else {
        pval <- 2 * pt(-abs(tstat), df=x$df.residual) #df here is ad-hoc
    }
    fixed <- cbind(est, est.se, tstat, pval)
    colnames(fixed) <- c("Estimate", "Std. Error", "t value", "Pr(> |t|)")
    rownames(fixed) <- names

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

    cat("\nRandom effects:\n")
    print(random, digits = digits, na.print = "NA", ...)

    cat("\nFixed effects:\n")
    printCoefmat(fixed, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)

    cat("\n(Dispersion parameter for ", x$family$family,
        " family taken to be ", format(x$dispersion, digits=digits), ")\n",
        sep="")

    cat("\n")
}
