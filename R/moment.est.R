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



# Compute 'sufficient statitics' for estimating fixed coef mean using moment method.
# Usually it will be used separately on each cores.
moment.est.mean.mapper <- function(coefficients, nfixed, subspace, precision, dispersion, start.cov = NULL)
{


  ngroups <- nrow(coefficients)
  dim <- ncol(coefficients)
  nrandom <- dim - nfixed
  names <- colnames(coefficients)

  fixed <- seq_len(nfixed)
  random  <- nfixed + seq_len(nrandom)

  if (is.null(start.cov))
    start.cov <- diag(nrandom)

  # compute mean estimate and covariance bias correction
  weight11.sum <- matrix(0,  nfixed, nfixed)
  weight1.coef.sum <- matrix(0,  nfixed, 1)

  for (i in seq_len(ngroups)) {
    u <- subspace[[i]]
    u1 <- u[fixed,,drop=FALSE]
    u2 <- u[random,,drop=FALSE]
    l <- precision[[i]]

    sigma2 <- dispersion[i]
    r <- length(l)

    if (r == 0L) {
      next
    }

    s <- sqrt(l)
    us <- u %*% diag(s, r, r)
    u1s <- us[fixed,,drop=FALSE]
    u2s <- us[random,,drop=FALSE]

    cov22 <- t(u2s) %*% start.cov %*% u2s
    w.inv <- cov22 + diag(sigma2, r, r)
    R.w.inv <- chol(w.inv)
    R.u1s.t <- backsolve(R.w.inv, t(u1s), transpose=TRUE)
    R.u2s.t <- backsolve(R.w.inv, t(u2s), transpose=TRUE)

    w11 <- t(R.u1s.t) %*% R.u1s.t
    w12 <- t(R.u1s.t) %*% R.u2s.t

    R2.u2s.t <- backsolve(R.w.inv, R.u2s.t)

    w1b <- (w11 %*% coefficients[i,fixed]
	    + w12 %*% coefficients[i,random])

    weight11.sum <- weight11.sum + w11
    weight1.coef.sum <- weight1.coef.sum + w1b
  }

  return(list(ngroups = ngroups, 
	      weight11.sum  = weight11.sum ,
	      weight1.coef.sum  = weight1.coef.sum ))

}

# Estimate fixed effect's coefficient. 
# Usually used as the 'combine' function in foreach function.
moment.est.mean.reducer <- function(...)
{

  ret <- list(...)
  if(length(ret)==0)
    stop("The input is empty.")


  ngroups <- 0
  wtot <- 0
  meanSum <- 0

  for(i in seq_len(length(ret))){
    ngroups <- ngroups + ret[[i]]$ngroups
    wtot <- wtot + ret[[i]]$weight11.sum
    meanSum <- meanSum + ret[[i]]$weight1.coef.sum
  }
  wtot <- wtot / ngroups
  meanSum <- meanSum / ngroups

  mean <- pseudo.solve(wtot, meanSum)
  if (attr(mean, "deficient")) {
    warning("cannot solve fixed effect moment equation due to rank deficiency")
  }
  mean.cov <- pseudo.solve(wtot) / ngroups
  attr(mean, "deficient") <- attr(mean.cov, "deficient") <- NULL

  est.mean <- list(mean = mean, mean.cov = mean.cov)
  return(est.mean)
}

# Compute 'sufficient statitics' for estimating ranef cov using moment method.
# Usually it will be used separately on each cores.
moment.est.cov.mapper <- function(coefficients, nfixed, subspace, precision, dispersion,
                       start.cov = NULL, diagcov = FALSE, mean)
{

  dim <- ncol(coefficients)
  ngroups <- nrow(coefficients)
  nrandom <- dim - nfixed
  fixed <- seq_len(nfixed)
  random  <- nfixed + seq_len(nrandom)

  names <- colnames(coefficients)

  if (is.null(start.cov))
    start.cov <- diag(nrandom)

  # compute mean estimate and covariance bias correction
  if(diagcov){
    bias.sum <- matrix(0, nrandom, 1)
    weight22.2.sum <- matrix(0,  nrandom, nrandom)
    weight22.coef.2.sum <- matrix(0, nrandom, 1)
  } else {
    bias.sum <- matrix(0,  nrandom, nrandom)
    weight22.2.sum <- matrix(0,  nrandom^2, nrandom^2)
    weight22.coef.2.sum <- matrix(0,  nrandom, nrandom)
  }

  for (i in seq_len(ngroups)) {
    u <- subspace[[i]]
    u1 <- u[fixed,,drop=FALSE]
    u2 <- u[random,,drop=FALSE]
    l <- precision[[i]]

    sigma2 <- dispersion[i]
    r <- length(l)
  
    if (r == 0L) {
      next
    }

    s <- sqrt(l)
    us <- u %*% diag(s, r, r)
    u1s <- us[fixed,,drop=FALSE]
    u2s <- us[random,,drop=FALSE]

    cov22 <- t(u2s) %*% start.cov %*% u2s
    w.inv <- cov22 + diag(sigma2, r, r)
    R.w.inv <- chol(w.inv)
    R.u1s.t <- backsolve(R.w.inv, t(u1s), transpose=TRUE)
    R.u2s.t <- backsolve(R.w.inv, t(u2s), transpose=TRUE)

    w12 <- t(R.u1s.t) %*% R.u2s.t
    w22 <- t(R.u2s.t) %*% R.u2s.t

    R2.u2s.t <- backsolve(R.w.inv, R.u2s.t)

    if(diagcov)
      B <- sigma2 * apply((R2.u2s.t)^2,2,sum)
    else 
      B <- sigma2 * t(R2.u2s.t) %*% R2.u2s.t

    bias.sum <- bias.sum + B



    if(diagcov){
      diff <- (t(w12) %*% (coefficients[i,fixed] - mean)
	       + w22 %*% coefficients[i,random])

      weight22.2.sum <- weight22.2.sum + w22^2
      weight22.coef.2.sum  <- weight22.coef.2.sum + diff^2
    } else {

      weight22.2.sum <- weight22.2.sum + kronecker(w22, w22)

      diff <- (t(w12) %*% (coefficients[i,fixed] - mean)
	       + w22 %*% coefficients[i,random])
      weight22.coef.2.sum <- weight22.coef.2.sum + diff %*% t(diff)
    }

  }

  ret <- list(ngroups = ngroups,
	      weight22.2.sum = weight22.2.sum, 
	      weight22.coef.2.sum = weight22.coef.2.sum,
	      bias.sum = bias.sum)
  return(ret)

}


# Estimate ranef covariance matrix, if control$diagcov is TRUE, i.e. only estimate diagonal element in covariance matrix.
# Usually used as the 'combine' function in foreach function.
moment.est.cov.reducer.diag<- function(...)
{

  ret <- list(...)
  if(length(ret)==0)
    stop("The input is empty.")

  ngroups <- 0
  wtot2 <- 0
  wt.cov <- 0
  wt.bias <- 0

  for(i in seq_len(length(ret))){
    ngroups <- ngroups + ret[[i]]$ngroups
    wtot2 <- wtot2 + ret[[i]]$weight22.2.sum
    wt.cov <- wt.cov + ret[[i]]$weight22.coef.2.sum
    wt.bias <- wt.bias + ret[[i]]$bias.sum
  }

  nrandom <- nrow(wt.bias)

  LHSinv <- solve(wtot2)
  diag_1 <- LHSinv %*% wt.cov
  diag_2 <- LHSinv %*% wt.bias
#  gamma <- min(min(diag_1 / diag_2),1)
#  cov <- diag( c(diag_1 - gamma * diag_2) , nrow = nrandom)
  cov <- diag( diag_1 -  diag_2 , nrow = nrandom)

  cov <- proj.psd(cov)  # ensure positive definite
  if (attr(cov, "modified") )
    warning(paste("moment-based covariance matrix estimate is not positive"
                  , " semi-definite; using projection"
                  , sep=""))
  attr(cov, "modified") <- NULL
  cov <- matrix(cov,nrow = nrandom)

  return(cov)

}

# Estimate ranef covariance matrix, if control$diagcov is FALSE, i.e. estimate every element in covariance matrix. 
# Usually used as the 'combine' function in foreach function.
moment.est.cov.reducer.exact <- function(...)
{

  ret <- list(...)
  if(length(ret)==0)
    stop("The input is empty.")


  ngroups <- 0
  wtot2 <- 0
  wt.cov <- 0
  wt.bias <- 0

  for(i in seq_len(length(ret))){
    ngroups <- ngroups + ret[[i]]$ngroups
    wtot2 <- wtot2 + ret[[i]]$weight22.2.sum
    wt.cov <- wt.cov + ret[[i]]$weight22.coef.2.sum
    wt.bias <- wt.bias + ret[[i]]$bias.sum
  }

  wtot2 <- wtot2/ngroups
  wt.cov <- wt.cov/ngroups
  wt.bias <- wt.bias/ngroups

  nrandom <- nrow(wt.bias)

  # construct an orthonormal basis for the space of symmetric
  # matrices
  q <- nrandom
  F <- matrix(0, q^2, q * (q + 1) / 2)
  j <- 0
  for (k in seq_len(q)) {
    for (l in seq_len(k)) {
      j <- j + 1
      f <- matrix(0, q, q)
      if (k == l) {
	f[k,l] <- 1
      } else {
	f[k,l] <- 1/sqrt(2)
	f[l,k] <- 1/sqrt(2)
      }
      F[,j] <- as.vector(f)
    }
  }

  # solve the moment equations
  tF.wtot2.F <- t(F) %*% wtot2 %*% F
  cov.vec <- pseudo.solve(tF.wtot2.F, t(F) %*% as.vector(wt.cov))
  bias.vec <- pseudo.solve(tF.wtot2.F, t(F) %*% as.vector(wt.bias))

  if (attr(cov.vec, "deficient")) {
    warning("cannot solve covariance moment equation due to rank deficiency")
  }

  # change back to original space
  cov <- matrix(F %*% cov.vec, nrandom, nrandom)
  bias <- matrix(F %*% bias.vec, nrandom, nrandom)

  # remove asymmetry arising from numerical errors
  cov <- 0.5 * (cov + t(cov))
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
  cov.adj <- (t(u.bias.t)
	      %*% diag((scale - l.bias) / scale, length(l.bias))
	      %*% u.bias.t)

  cov <- proj.psd(cov.adj)  # ensure positive definite
  if (attr(cov, "modified") || length(l) < nrow(cov) || scale != 1)
    warning(paste("moment-based covariance matrix estimate is not positive"
		  , " semi-definite; using projection"
		  , sep=""))
  attr(cov, "modified") <- NULL
  return(cov)
}


moment.est <- function(coefficients, nfixed, subspace, precision, dispersion,
                       start.cov = NULL, parallel = FALSE, diagcov = FALSE)
{

  logging::loginfo("Estimating moments", logger="mbest.mhglm.fit")
  ngroups <- nrow(coefficients)
  dim <- ncol(coefficients)
  names <- colnames(coefficients)
  if (ngroups == 0L || dim == 0L){
    cov <- matrix(0, dim, dim)
    dimnames(cov) <- list(names, names)
    return(cov)
  }


  logging::loginfo("Computing mean estimate", logger="mbest.mhglm.fit")
  i <- NULL
  # fixef
  if(parallel){
    est.mean <- foreach(i=seq_len(ngroups),
			.combine = 'moment.est.mean.reducer',
			.multicombine=TRUE ) %dopar% {
      moment.est.mean.mapper(coefficients[i,,drop = FALSE], nfixed, 
			     list(subspace[[i]]),list(precision[[i]]),
			     dispersion[i], start.cov = start.cov) 
    }
  } else {
    mean.info <- moment.est.mean.mapper(coefficients, nfixed, subspace,precision,
					dispersion, start.cov = start.cov) 
    est.mean <- moment.est.mean.reducer(mean.info)
  }

  logging::loginfo("Computing covariance estimate", logger="mbest.mhglm.fit")
  # ranef
  if(parallel){
    if(diagcov){
      est.cov <- foreach(i=seq_len(ngroups), .combine = 'moment.est.cov.reducer.diag', .multicombine = TRUE) %dopar%{
	moment.est.cov.mapper(coefficients[i,,drop = FALSE],nfixed,
			      list(subspace[[i]]),list(precision[[i]]),dispersion[i],
			      start.cov = start.cov, diagcov = diagcov, 
			      est.mean$mean) }
    } else {
      est.cov <- foreach(i=seq_len(ngroups), .combine = 'moment.est.cov.reducer.exact', .multicombine = TRUE) %dopar%{
	moment.est.cov.mapper(coefficients[i,,drop = FALSE], nfixed, 
			      list(subspace[[i]]),list(precision[[i]]),dispersion[i],
			      start.cov = start.cov, diagcov = diagcov, 
			      est.mean$mean) }
    }
  } else {
    cov.info <- moment.est.cov.mapper(coefficients, nfixed,
				      subspace,precision,dispersion,
				      start.cov = start.cov, diagcov = diagcov, 
				      est.mean$mean) 
    if(diagcov){
      est.cov <- moment.est.cov.reducer.diag(cov.info)
    } else {
      est.cov <- moment.est.cov.reducer.exact(cov.info)
    }
  }

  list(mean=est.mean$mean, mean.cov=est.mean$mean.cov, cov=est.cov)

}


