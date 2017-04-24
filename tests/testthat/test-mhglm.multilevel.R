context("mhglm_ml")
library("lme4")

expect_equal_tol4 <- function(...) expect_equal(..., tolerance = 1e-4)

test_that("same results for mhglm and mhglm_ml on sleepstudy", {
    model <- mhglm(Reaction ~ Days + (Days | Subject), data = sleepstudy)
    model_ml <- mhglm_ml(Reaction ~ Days + (Days | Subject), data = sleepstudy)

    expect_equal(family(model_ml), family(model))
    expect_equal(fitted(model_ml), fitted(model))
    expect_equal(fixef(model_ml), fixef(model))
    expect_equal(model.frame(model_ml), model.frame(model))
    expect_equal(model.matrix(model_ml), model.matrix(model))
    expect_equal(predict(model_ml), predict(model))
    expect_equal(residuals(model_ml), residuals(model))
    expect_equal(sigma(model_ml), sigma(model))
    expect_equal(terms(model_ml), terms(model))
    expect_equal(VarCorr(model_ml), VarCorr(model))
    expect_equal(vcov(model_ml), vcov(model))
    expect_equal(weights(model_ml), weights(model))

    expect_equivalent(ranef(model_ml), ranef(model)) # class is different
})

test_that("same results for mhglm and mhglm_ml on sleepstudy, in parallel", {
    model <- mhglm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                   control = list(parallel = TRUE))
    model_ml <- mhglm_ml(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                         control = list(parallel = TRUE))

    expect_equal(family(model_ml), family(model))
    expect_equal(fitted(model_ml), fitted(model))
    expect_equal(fixef(model_ml), fixef(model))
    expect_equal(model.frame(model_ml), model.frame(model))
    expect_equal(model.matrix(model_ml), model.matrix(model))
    expect_equal(predict(model_ml), predict(model))
    expect_equal(residuals(model_ml), residuals(model))
    expect_equal(sigma(model_ml), sigma(model))
    expect_equal(terms(model_ml), terms(model))
    expect_equal(VarCorr(model_ml), VarCorr(model))
    expect_equal(vcov(model_ml), vcov(model))
    expect_equal(weights(model_ml), weights(model))

    expect_equivalent(ranef(model_ml), ranef(model)) # class is different
})

test_that("succeeds on simdata for two levels", {
    model <- mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g2), data = simdata)

    # fixef
    fixef0 <- c("(Intercept)" = 0.8600, "x" = 0.4713)
    fixef_m <- fixef(model)
    expect_equal_tol4(fixef_m, fixef0)

    # vcov
    vcov0 <- matrix(c(0.1812, 0.0252, 0.0252, 0.0798), nrow = 2L, ncol = 2L,
        dimnames = list(c("(Intercept)", "x"), c("(Intercept)", "x")))
    vcov_m <- vcov(model)
    expect_equal_tol4(vcov_m, vcov0)

    # VarCorr
    varcor0 <- list(
        g1 = matrix(c(1.7208, 0.2431, 0.2431, 0.6780), nrow = 2L, ncol = 2L,
            dimnames = list(c("(Intercept)", "x"), c("(Intercept)", "x"))),
        g2 = matrix(c(0.3336, 0.0014, 0.0014, 0.5170), nrow = 2L, ncol = 2L,
            dimnames = list(c("(Intercept)", "x"), c("(Intercept)", "x"))),
        sc = 1.0343)
    varcor_m <- VarCorr(model)

    expect_equal(attr(varcor_m$g1, "stddev"), sqrt(diag(varcor_m$g1)))
    expect_equivalent(attr(varcor_m$g1, "correlation"),
                      varcor_m$g1 / tcrossprod(attr(varcor_m$g1, "stddev")))
    expect_equivalent(round(varcor_m$g1, 4), varcor0$g1)

    expect_equal(attr(varcor_m$g2, "stddev"), sqrt(diag(varcor_m$g2)))
    expect_equivalent(attr(varcor_m$g2, "correlation"),
                      varcor_m$g2 / tcrossprod(attr(varcor_m$g2, "stddev")))
    expect_equivalent(round(varcor_m$g2, 4), varcor0$g2)

    expect_equal_tol4(attr(varcor_m, "sc"), varcor0$sc)

    ranef0 <- list(
        g1 = structure(list(
            `(Intercept)` = c(0.5759, 0.2531, 0.1683, -1.3702, 0.2443, -1.8005,
                0.5905, -1.5495, 0.0243, 2.8659),
            x = c(0.216, -1.6434, 1.2734, -0.1194, -0.5528, -0.3414, 0.9547,
                -0.2058, 0.0712, 0.339)),
            row.names = as.character(1:10),
            .Names = c("(Intercept)", "x"),
            class = "data.frame"),
        g2 = structure(list(
            `(Intercept)` =
                c(-0.1039, -0.5185, 0.3325, 0.0331, 0.3228, -0.4613,
                0.2889, -0.1203, -0.0951, 0.6658, -0.0785, -0.6297, 0.4605,
                1.1001, -1.3396, 0.5025, 0.2419, 0.4263, -0.4063, -1.1126,
                -0.7274, -0.8474, -0.0537, 1.0752, 0.424, 0.7172, -0.4676,
                0.7311, -1.1279, -0.4519, -0.2108, 0.6348, 0.641, -0.4896,
                -0.4572, 0.5019, -0.314, 0.1451, -0.3624, -0.0813, 0.285,
                -0.6433, -0.1163, -0.3561, 0.818, 0.597, 0.0872, -0.0127,
                -0.5505, 0.1268),
            x = c(-1.0476, -0.4643, 0.1117, 1.1567, 0.39, -0.1245, -0.5155,
                0.7865, -0.5186, -0.3329, 1.5228, -0.4958, -0.5501, 0.0928,
                0.762, -0.478, -0.6605, 0.0643, 0.5636, 0.7897, 0.1156, -0.6142,
                -1.2659, 0.5212, 0.7736, 0.2684, -0.2731, 0.5588, -0.358,
                -0.5415, 1.2585, 2.0112, -1.0565, 0.9565, -0.913, -0.4056,
                0.209, 0.015, 0.4263, -0.1617, 0.4982, -0.0212, 0.3323, -0.3245,
                -0.5071, 0.6005, -0.7002, -0.2273, 0.89, -0.7108)),
            row.names = as.character(1:50),
            .Names = c("(Intercept)", "x"),
            class = "data.frame"))
    ranef_m <- ranef(model)

    expect_equal_tol4(ranef_m$g1, ranef0$g1)
    expect_equal_tol4(ranef_m$g2, ranef0$g2)
})

test_that("succeeds on multiple 2-level formula permutations", {
    model1 <- mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g2), data = simdata)
    model2 <- mhglm_ml(y ~ 1 + x + (1 + x | g2) + (1 + x | g1), data = simdata)
    model3 <- mhglm_ml(y ~ 1 + x + (1 + x | g1/g2), data = simdata)
    model4 <- mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g1:g2),
                       data = simdata)

    expect_equal(fixef(model1), fixef(model2))
    expect_equal(fixef(model1), fixef(model3))
    expect_equal(fixef(model1), fixef(model4))

    expect_equivalent(ranef(model1), ranef(model2))
    expect_equivalent(ranef(model1), ranef(model3))
    expect_equivalent(ranef(model1), ranef(model4))
})

test_that("succeeds on multiple 3-level formula permutations", {
    model1 <- mhglm_ml(y ~ x + (x | g1) + (x | g2) + (x | g3), data = simdata)
    model2 <- mhglm_ml(y ~ x + (x | g1/g2/g3), data = simdata)
    model3 <- mhglm_ml(y ~ x + (x | g1/g2) + (x | g3), data = simdata)
    model4 <- mhglm_ml(y ~ x + (x | g1/g2) + (x | g1:g2:g3), data = simdata)
    model5 <- mhglm_ml(y ~ x + (x | g1) + (x | g2/g3), data = simdata)
    model6 <- mhglm_ml(y ~ x + (x | g1) + (x | g1:g2) + (x | g1:g2:g3),
                       data = simdata)
    model7 <- mhglm_ml(y ~ x + (x | g1:g2:g3) + (x | g1:g2) + (x | g1),
                       data = simdata)

    expect_equal(fixef(model1), fixef(model2))
    expect_equal(fixef(model1), fixef(model3))
    expect_equal(fixef(model1), fixef(model4))
    expect_equal(fixef(model1), fixef(model5))
    expect_equal(fixef(model1), fixef(model6))
    expect_equal(fixef(model1), fixef(model7))

    expect_equivalent(ranef(model1), ranef(model2))
    expect_equivalent(ranef(model1), ranef(model3))
    expect_equivalent(ranef(model1), ranef(model4))
    expect_equivalent(ranef(model1), ranef(model5))
    expect_equivalent(ranef(model1), ranef(model6))
    expect_equivalent(ranef(model1), ranef(model7))
})

test_that("success using diagonal covariance", {
    model <- mhglm_ml(y ~ x + (x | g1) + (x | g2) + (x | g3), data = simdata,
                      control = mhglm_ml.control(diagcov = TRUE))
    varcor <- VarCorr(model)
    expect_equal(varcor$g1[1, 2], 0)
    expect_equal(varcor$g2[1, 2], 0)
    expect_equal(varcor$g3[1, 2], 0)
})
