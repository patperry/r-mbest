# 1.glmfit.R
# (Map step) For each piece of data indexed by i
# - fit glm to it
# - compute dispersion.sum, df.residual.tot, coef, Rp, subspace, precision
# - write out info.i.rds, subspace.i.rds, disp.i.csv
#
# (Reduce step) Pool information together
# - read in all disp.*.csv
# - compute dispersion.tot, df.residual.tot


start <- etastart <- mustart <- intercept <- TRUE
control <- do.call("mhglm.control", control)
nfiles <- length(datafile_list)

for (dataid in seq_len(nfiles)) {
    # load one piece of data
    datafile <- paste0(datafile_folder, datafile_list[dataid])
    data <- readRDS(datafile)
    x <- data$x
    y <- data$y
    z <- x
    group <- data$group

    # set parameters
    nobs <- NROW(y)
    weights <- rep(1, nobs)
    offset <- rep(0, nobs)

    # save data specific information
    infofilename <- paste0(datafile_folder, "info.", dataid, ".rds")
    saveRDS(list(nobs = nobs, weights = weights, offset = offset, y = y,
                 group = group), file = infofilename)

    # start mhglm.fit
    xnames <- dimnames(x)[[2L]]
    znames <- dimnames(z)[[2L]]
    ynames <- if (is.matrix(y)) rownames(y) else names(y)
    nfixed <- ncol(x)
    nrandom <- ncol(z)
    nvars <- nfixed + nrandom
    fixed <- (seq_len(nvars) <= nfixed)
    random <- !fixed

    # doesn't support standardized
    R.fixed <- diag(nfixed)
    attr(R.fixed, "pivot") <- seq_len(nfixed)
    attr(R.fixed, "rank") <- nfixed
    R.random <- diag(nrandom)
    attr(R.random, "pivot") <- seq_len(nrandom)
    attr(R.random, "rank") <- nrandom

    rank.fixed <- attr(R.fixed, "rank")
    rank.random <- attr(R.random, "rank")
    rank <- rank.fixed + rank.random
    r1.fixed <- seq_len(rank.fixed)
    r1.random <- seq_len(rank.random)
    r1 <- seq_len(rank)
    pivot.fixed <- which(fixed)[attr(R.fixed, "pivot")]
    pivot.random <- which(random)[attr(R.random, "pivot")]
    pivot <- c(pivot.fixed[r1.fixed], pivot.random[r1.random],
    pivot.fixed[-r1.fixed], pivot.random[-r1.random])
    R <- matrix(0, rank, rank)
    R[r1.fixed, r1.fixed] <- R.fixed[r1.fixed, r1.fixed]
    (R[rank.fixed + r1.random, rank.fixed + r1.random] <-
            R.random[r1.random, r1.random])


    if (!is.null(start)) {
        if (length(start) != nfixed) {
            stop(gettextf(
                paste("length of 'start' should equal %d and correspond to",
                      "initial coefs for %s"),
                nfixed, paste(deparse(xnames), collapse = ", ")), domain = NA)
        }
        start <- c(start, numeric(nrandom))
    }


    # fit group specific glm
    m <- rdglm.group.fit(x = cbind(x, z), y = y, group = group,
                         weights = weights, start = start, etastart = etastart,
                         mustart = mustart, offset = offset, family = family,
                         parallel = control$parallel,
                         control = control$fit.control,
                         method = control$fit.method, intercept = intercept)

    # compute and save dispersion information
    dispfilename <- paste0(datafile_folder, "disp.", dataid, ".csv")

    df.residual.tot <- sum(m$df.residual)
    dispersion.sum <- sum(m$df.residual * m$dispersion)
    ngroups <- m$ngroups

    estdispersion <- data.frame(ngroups, dispersion.sum, df.residual.tot)
    write.table(estdispersion, file = dispfilename, sep = ",", quote = FALSE,
                col.names = FALSE, row.names = FALSE)

    # compute group-sqecific precision square root (unpivoted)
    Rp <- as.list(rep(NULL, m$ngroups))
    if (control$parallel) {
        Rp <- foreach(i = seq_len(m$ngroups)) %dopar% {
            qr.i <- m$qr[[i]]
            return(unpivotRp(qr.i, nvars))
        }
    } else {
        if (i %in% seq_len(m$ngroups)) {
            qr.i <- m$qr[[i]]
            Rp[[i]] <- unpivotRp(qr.i, nvars)
        }
    }
    names(Rp) <- names(m$qr)

    # compute standardized coeficients:
    coef <- m$coefficients[, pivot[r1], drop = FALSE] %*% t(R)

    # compute group-specific subspace and precisions
    if (control$parallel) {
        results <- foreach(i = seq_len(m$ngroups)) %dopar% {
            if (nrow(Rp[[i]]) > 0L) {
                prec.sqrt <- backsolve(R, t(Rp[[i]][, pivot[r1], drop = FALSE]),
                transpose = TRUE)
                prec.sqrt.svd <- svd(prec.sqrt)
                return(list(subspace = prec.sqrt.svd$u,
                            precision = (prec.sqrt.svd$d) ^ 2))
            } else {
                return(list(subspace = matrix(0, rank, 0),
                            precision = numeric()))
            }
        }
        subspace <- lapply(results, function(x) x[["subspace"]])
        precision <- lapply(results, function(x) x[["precision"]])
        rm(results)
    } else {
        subspace <- list()
        precision <- list()
        for (i in seq_len(m$ngroups)) {
            if (nrow(Rp[[i]]) > 0L) {
                prec.sqrt <- backsolve(R, t(Rp[[i]][, pivot[r1], drop = FALSE]),
                transpose = TRUE)
                prec.sqrt.svd <- svd(prec.sqrt)
                subspace[[i]] <- prec.sqrt.svd$u
                precision[[i]] <- (prec.sqrt.svd$d) ^ 2
            } else {
                subspace[[i]] <- matrix(0, rank, 0)
                precision[[i]] <- numeric()
            }
        }
    }
    names(subspace) <- names(Rp)
    names(precision) <- names(Rp)

    # save group specific coef, Rp, subspace, precision
    subspacefilename <- paste0(datafile_folder, "subspace.", dataid, ".rds")
    saveRDS(list(ngroups = m$ngroups, coef = coef, Rp = Rp,
    subspace = subspace, precision = precision), file = subspacefilename)
}


#----------
# pool dispersion information togeter and estimate global dispersion
dispfiles <- grep("disp", list.files(datafile_folder, full.names = TRUE),
                  value = TRUE)
disp.summary <- data.frame()
for (i in seq_len(length(dispfiles))) {
  disp.summary <- rbind(disp.summary,
			read.table(dispfiles[i],sep = ","))
}

disp.summary <- colSums(disp.summary)
names(disp.summary) <- c("ngroups", "dispersion.sum", "df.residual.tot")

dispersion.tot <- if (disp.summary["df.residual.tot"] > 0) {
    disp.summary["dispersion.sum"] / disp.summary["df.residual.tot"]
} else {
    disp.summary["dispersion.sum"] / disp.summary["ngroups"]
}

df.residual.tot <- disp.summary["df.residual.tot"]
