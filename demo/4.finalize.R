# 6.finalize.R
# Construct the same return object as mhglm.fit.
# Save the fit object to mhglmfit.rds.

#--------------------
# combind subspace, precision, 
# y, group, weights, offset, nobs 
subspace <- c()
precision <- c()
y <- c()
group <- c()
weights <- c()
offset <- c()
nobs <- 0
for(dataid in seq_len(nfiles)){

  infofilename <- paste0(datafile_folder,'info.',dataid,".rds")
  info <- readRDS(infofilename)
  subspacefile <- paste0(datafile_folder,'subspace.',dataid,'.rds')
  s <- readRDS(subspacefile)

  subspace <- c(subspace, s$subspace)
  precision <- c(precision, s$precision)
  y <- c(y,info$y)
  group <- c(group,list(info$group))
  weights <- c(weights, info$weights)
  offset <- c(offset, info$offset)
  nobs <- nobs + info$nobs
}
group <- unlist(group)


#--------------------
# change back to original coordinates
mean <- est$mean
mean.cov <- est$mean.cov
cov <- est$cov

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

#--------------------
# set coordinate names
names(coef.mean) <- xnames
dimnames(coef.mean.cov) <- list(xnames, xnames)
dimnames(coef.cov) <- list(znames, znames)

#--------------------
# construct the fit object and save to disk
fit <- list(family = family, coefficient.mean = coef.mean,
	    coefficient.mean.cov = coef.mean.cov,
	    coefficient.cov = coef.cov, coefficients = coef,
	    subspace = subspace, precision = precision,
	    dispersion = dispersion.tot, df.residual = df.residual.tot,
	    R = R, rank = rank, rank.fixed = rank.fixed,
	    rank.random = rank.random, pivot = pivot, y = y, group = group,
	    prior.weights = weights, offset = offset, nobs = nobs)


fitfile <- paste0(datafile_folder,"mhglmfit.rds")
saveRDS(fit,file = fitfile)

