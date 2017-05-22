# 3.fixef.R 
# Read in momentest.rds (contains estimate for ranef covariance), if it exists. 
# Otherwise set the initial covariance estimate to NULL.
# 
# (Map step) For each piece of data indexed by i
# - read in subspace.i.rds
# - compute fixef's sufficient statistics
# - write out fixef.i.rds
# 
# (Reduce step) Pool information together
# - read in all fixef.i.rds 
# - estimate fixef

#---------- 
# read in initial covariance estimate.
# if no initial estimate is provided, set to NULL.
mestfile<- paste0(datafile_folder,'momentest.rds',sep = "")
if(file.exists(mestfile)){
  est0 <- readRDS(mestfile)
} else {
  est0 <- NULL
} 
start.cov <- est0$cov

#---------- 
# compute and save sufficient statistics for estimating fixef
for(dataid in seq_len(nfiles)){

  subspacefile <- paste0(datafile_folder,'subspace.',dataid,'.rds')
  s <- readRDS(subspacefile)

  dispersion <- rep(dispersion.tot, s$ngroups)
  fixef.sep <- moment.est.mean.mapper(s$coef, rank.fixed, s$subspace, s$precision, dispersion, start.cov)

  fixeffilename <- paste(datafile_folder,'fixef.',dataid,'.rds',sep = "")
  saveRDS(fixef.sep, fixeffilename)
}

#---------- 
# estimate fixef
fixeffiles <- grep('fixef.*.rds',list.files(datafile_folder,full.names = TRUE),value = TRUE)
mean.info <- as.list(rep(NULL,nfiles))

for (i in seq_len(nfiles)){
  fixef.ret <- readRDS(fixeffiles[i])
  mean.info[[i]] <- fixef.ret
}

est.mean <- moment.est.mean.reducer(mean.info)


