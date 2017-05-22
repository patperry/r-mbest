library(devtools)
load_all("../R/", export_all = TRUE) # export all functions, for convenience

#----------
# generate simulation samples;
# split data by group and write to separate files
#----------
family <- gaussian()
datafile_folder <- "./"
source("0.sim.R")


#----------
# The following lines of code demonstrate how to run
# mhglm.fit on very large dataset.
#
# The main idea is to process one trunk of data at
# a time, save the necessary information and move on to
# the next trunk. After one pass of data, combine the
# saved information to make final estimate. This is the
# classic Map/Reduce method.
#
# This method requires that data from same group must be
# processed together, therefore should be saved in same file.
#
# The detailed procedures are as below.
# 1. Set parameters. In addition to parameters like control and family,
# one also need to specify where the data lives and the filename list.
#
# 2. Fit glm to each trunk of data. For each trunk of data, run rdglm.group.fit,
# and then compute dispersion.sum, df.residual.tot, Rp, coef,
# subspace, precision. Save the above information for each trunk.
# Once all the data are processed, combine dispersion.sum and df.residual
# to estimate the overall dispersion.
#
# 3. Compute coefficient mean and covariance estimates.
# - 2.fixed.R: for each trunk of data, compute and save the sufficient statistics for
# computing fixef (coefficient mean). Once all the data are processed, combine the
# saved statistics to estimate fixef.
# - 3. ranef.R: for each trunk of data, compute and save the sufficient statistics for
# computing ranef covariance matrix. Once all the data are processed, combine the
# saved statistics to estimate ranef.cov.
# - iterate over the above two steps.
#
# 4. Save the same output as mhglm.fit. Construct the same object as mhglm.fit would return.
# Save to disk.
#----------

# 1. Set parameters
control <- list(diagcov = TRUE, standardize = FALSE)
family <- gaussian()

# where are the partitioned data files
datafile_folder <- "./"

# the list of data files to process
datafile_list <- grep("data", list.files(datafile_folder), value = TRUE)


# 2. Fit glm to each trunk of data.
source("1.glmfit.R")

# 3. Compute coefficient mean and covariance estimates
iter <- 0
while (iter < control$steps) {
  source("2.fixef.R")
  source("3.ranef.R")
  iter <- iter + 1
}
source("2.fixef.R")
source("3.ranef.R")

# 4. Save the same output as mhglm.fit
source("4.finalize.R")
