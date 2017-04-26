library(mbest)

parse_formula <- function
## Parse the formula and return a list of objects that
## will be fed into fitting function.
##
## Input: formula and data (N samples).
## Output:
## - x: fixed effects matrix of N x p (p is dimension of fixed effects)
## - y: vector of response, length N 
## - z: list of d components (d is number of grouping levels)
##  - z[[1]]: first level random effects matrix of N x q1 (q1 is dimension of random effects on first level)
##  - z[[2]]: second level ...
##  - etc.
## - group: list of d components (d is number of grouping levels)
##  - group[[1]]: factor vector of length N. Each factor corresoponds to a group on first level
##  - group[[2]]: factor vector of length N. Each factor corresoponds to a (nested) group on second level. 
##  - etc.
(formula, 
 data, 
 method = NULL)
{
    # model frame
    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf$formula <- lme4::subbars(eval(mf$formula))
    mf <- eval(mf, parent.frame())

    # method
    if (identical(method, "model.frame"))
        return(mf)

    # terms
    mt <- attr(mf, "terms")

    # response
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }

    mt.fixed <- delete.response(terms(lme4::nobars(formula), data=data))
    X <- if (!is.empty.model(mt.fixed)){
        model.matrix(mt.fixed, mf, contrasts)
    } else {matrix(, NROW(Y), 0L)}


    ## use customized findbars
    # so (1|a/b) -> (1|a) + (1|b:a)
    # lme::findbars: (1|a/b) -> (1|b:a) + (1|a)·
    bars <- findbars(formula)
    if (length(bars) >= 1L) {
      Z <- list()
      Group <- list()
      group.call <- list()
      mt.random <- list()

      for(i in seq_along(bars)){
          b <- bars[[i]]
          mf1 <- mf
          for (v in all.vars(b[[3L]])) {
              mf1[[v]] <- factor(mf1[[v]])
          }
          group.call.tmp <- substitute(factor(fac), list(fac = b[[3L]]))
          group.tmp <- eval(group.call.tmp, mf1)

          if (all(is.na(group.tmp)))
              stop("Invalid grouping factor specification, ", deparse(b[[3L]]))

          mt.random.tmp <- terms(eval(substitute(~trms, list(trms = b[[2L]]))),
                                 data = data)
          z.tmp <- if (!is.empty.model(mt.random.tmp)){
              model.matrix(mt.random.tmp, mf, contrasts)
          }else matrix(, NROW(Y), 0L)

          Z[[i]] <- z.tmp
          Group[[i]] <- group.tmp
          group.call[[i]] <- group.call.tmp
          mt.random[[i]] <- mt.random.tmp
      }
    } else { # length(bars) == 0L
        Group <- factor(character(NROW(Y)))
        mt.random <- terms(~ -1, data=data)
        Z <- matrix(, NROW(Y), 0L)
    }

    list(x = X, y = Y, z = Z, group = Group)

}

## Example
load('../data/simdata.rda')
# formula <- y ~ x + (x|g1/g2/g3)
formula <- y ~ x + (x|g1) + (x|g1:g2) + (x|g1:g2:g3)
data <- simdata
ret <- parse_formula(formula, data)

