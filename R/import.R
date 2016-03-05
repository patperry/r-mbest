## base packages
#' @importFrom stats .checkMFClasses, .getXlevels, coef, delete.response
#' @importFrom stats family, fitted, gaussian, is.empty.model, model.extract
#' @importFrom stats model.frame, model.offset, model.matrix, model.response
#' @importFrom stats model.weights, na.pass, naresid, nobs, pnorm, predict, 
#' @importFrom stats printCoefmat, pt, residuals, terms,update, vcov, weights
#' 

## Recommended packages
#' @importFrom lme4 findbars, fixef, nobars, ranef, subbars, VarCorr 
#' 

# Roxygen cannot handle conditional import!!!
#' @importFrom lme4, sigma

#' @useDynLib mbest
#' @importFrom Rcpp sourceCpp 