# 4.ranef.R
# (Map step) For each piece of data indexed by i
# - read in subspace.i.rds
# - compute ranef cov's sufficient statistics
# - write out ranef.i.rds
#
# (Reduce step) Pool information together
# - estimate ranef cov
# - write out momentest.rds

#----------
# compute and save sufficient statistics for estimating ranef cov
for (dataid in seq_len(nfiles)) {
    subspacefile <- paste0(datafile_folder, "subspace.", dataid, ".rds")
    s <- readRDS(subspacefile)

    dispersion <- rep(dispersion.tot, s$ngroups)
    ranef.sep <- moment.est.cov.mapper(s$coef, rank.fixed, s$subspace,
                                       s$precision, dispersion, start.cov,
                                       diagcov = control$diagcov,
                                       mean = est.mean$mean)

    raneffilename <- paste0(datafile_folder, "ranef.", dataid, ".rds")
    saveRDS(ranef.sep, raneffilename)
}

#----------
# estimate ranef cov
raneffiles <- grep("ranef.*.rds",
                   list.files(datafile_folder, full.names = TRUE),
                   value = TRUE)
cov.info <- list()
for (i in seq_len(nfiles)) {
    ranef.ret <- readRDS(raneffiles[i])
    cov.info[[i]] <- ranef.ret
}
est.cov <- moment.est.cov.reducer(cov.info, control$diagcov)

# save to disk
est <- list(mean = est.mean$mean, mean.cov = est.mean$mean.cov, cov = est.cov)
mestfile <- paste0(datafile_folder, "momentest.rds")
saveRDS(est, file = mestfile)
