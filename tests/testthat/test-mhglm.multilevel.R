context("mhglm_ml")
library("lme4")

g_data <- simdata(n = 30, m_per_level = c(10, 5, 2), sd_intercept = c(1, 1, 1),
                  sd_slope = c(1, 1, 1), family = "gaussian", seed = 12345)
b_data <- simdata(n = 30, m_per_level = c(10, 5, 2), sd_intercept = c(1, 1, 1),
                  sd_slope = c(1, 1, 1), family = "binomial", seed = 12345)

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
    model <- mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g2:g1), data = g_data)

    # fixef
    fixef0 <- c("(Intercept)" = 0.6947, "x" = 0.8033)
    fixef_m <- fixef(model)
    expect_equal_tol4(fixef_m, fixef0)

    # vcov
    vcov0 <- matrix(c(0.0635, 0.0135, 0.0135, 0.1314), nrow = 2L, ncol = 2L,
        dimnames = list(c("(Intercept)", "x"), c("(Intercept)", "x")))
    vcov_m <- vcov(model)
    expect_equal_tol4(vcov_m, vcov0)

    # VarCorr
    varcor0 <- list(
        g1 = matrix(c(0.4608, 0.1078, 0.1078, 1.0903), nrow = 2L, ncol = 2L,
            dimnames = list(c("(Intercept)", "x"), c("(Intercept)", "x"))),
        g2 = matrix(c(0.9467, 0.0032, 0.0032, 1.1192), nrow = 2L, ncol = 2L,
            dimnames = list(c("(Intercept)", "x"), c("(Intercept)", "x"))),
        sc = 1.3051)
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

    ranef0 <- structure(
        list(g1 = structure(
            list(`(Intercept)` = c(
                -0.0869, 0.5135, 0.2775, -0.9865, 0.6647, 0.0608, -0.3574,
                -0.7083, -0.3076, 0.933
            ),
            x = c(0.1599, -0.8793, -0.0248, -0.6779, 1.0913, -2.038, 1.168,
                  -0.2827, 0.7953, 0.6829)
            ),
            .Names = c("(Intercept)", "x"),
            row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
            class = "data.frame"
        ),
        `g2:g1` = structure(
            list(`(Intercept)` = c(
                1.1348, -0.6704, 0.3165, 0.2033, -1.3342, -0.8619, 1.0774,
                -0.8982, 0.4095, 1.1918, 0.7672, .0302, -0.9923, 2.5463,
                -1.1255, 0.3451, -0.428, -0.2005, -0.7753, 0.56, 1.3716, 1.1153,
                0.6209, 0.1745, -1.5804, -0.7036, 1.0323, -0.1678, 1.7321,
                -1.4476, -0.0917, 0.1559, 0.5558, -0.3264, -0.2809, 0.1841,
                0.8649, -1.5292, -1.0628, 0.271, -1.0486, 0.7977, -0.9432,
                0.7465, 0.2376, 1.3892, 2.6352, -0.7846, 0.6537, -0.5881
            ),
            x = c(
                -0.9315, 1.2815, -1.789, 0.4104, 1.4969, 0.4604, 0.9546,
                -1.0615, -0.161, -0.3434, -2.4954, 0.5899, 1.1695, -0.2297,
                0.3344, -0.2977, 0.2569, -0.7495, -0.0553, 0.6233, 0.2716,
                0.6725, 0.4516, -0.4217, -0.3452, -0.7267, 0.9873, 0.4423,
                -1.1059, -0.8408, 0.4285, 2.3327, 0.2993, 1.7577, -2.0216,
                0.9576, -0.4641, -1.5246, 1.0297, -0.3684, -0.0877, 1.4684,
                -0.8841, 0.112, -0.4388, 0.7136, -0.7475, -0.3823, 2.6744,
                -1.0387
            )),
            .Names = c("(Intercept)", "x"),
            row.names = c(
                "1:1", "2:1", "3:1", "4:1", "5:1", "6:2", "7:2", "8:2", "9:2",
                "10:2", "11:3", "12:3", "13:3", "14:3", "15:3", "16:4", "17:4",
                "18:4", "19:4", "20:4", "21:5", "22:5", "23:5", "24:5", "25:5",
                "26:6", "27:6", "28:6", "29:6", "30:6", "31:7", "32:7", "33:7",
                "34:7", "35:7", "36:8", "37:8", "38:8", "39:8", "40:8", "41:9",
                "42:9", "43:9", "44:9", "45:9", "46:10", "47:10", "48:10",
                "49:10", "50:10"
            ),
            class = "data.frame")),
        .Names = c("g1", "g2:g1"),
        class = "ranef.mhglm_ml")

    ranef_m <- ranef(model)

    expect_equal_tol4(ranef_m$g1, ranef0$g1)
    expect_equal_tol4(ranef_m$`g2:g1`, ranef0$`g2:g1`)
})

test_that("succeeds on multiple 2-level formula permutations", {
    model1 <- mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g2:g1),
                       data = g_data)
    model2 <- mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g1:g2),
                       data = g_data)
    model3 <- mhglm_ml(y ~ 1 + x + (1 + x | g2:g1) + (1 + x | g1),
                       data = g_data)
    model4 <- mhglm_ml(y ~ 1 + x + (1 + x | g1:g2) + (1 + x | g1),
                       data = g_data)
    model5 <- mhglm_ml(y ~ 1 + x + (1 + x | g1/g2), data = g_data)

    expect_equal(fixef(model1), fixef(model2))
    expect_equal(fixef(model1), fixef(model3))
    expect_equal(fixef(model1), fixef(model4))
    expect_equal(fixef(model1), fixef(model5))

    expect_equivalent(ranef(model1), ranef(model2))
    expect_equivalent(ranef(model1), ranef(model3))
    expect_equivalent(ranef(model1), ranef(model4))
    expect_equivalent(ranef(model1), ranef(model5))
})

test_that("fails on misspecified 2-level formulas", {
    expect_error(
        mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g2), data = g_data))
    expect_error(
        mhglm_ml(y ~ 1 + x + (1 + x | g2) + (1 + x | g1), data = g_data))
    expect_error(
        mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g1:g2) + (1 + x | g2),
                 data = g_data))
    expect_error(
        mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g1:g2) + (1 + x | g2:g1),
                 data = g_data))
})

test_that("succeeds on multiple 3-level formula permutations", {
    model1 <- mhglm_ml(y ~ x + (x | g1) + (x | g2:g1) + (x | g3:g2:g1),
                       data = g_data)
    model2 <- mhglm_ml(y ~ x + (x | g1/g2/g3), data = g_data)
    model3 <- mhglm_ml(y ~ x + (x | g1/g2) + (x | g3:g2:g1), data = g_data)
    model4 <- mhglm_ml(y ~ x + (x | g1/g2) + (x | g1:g2:g3), data = g_data)
    model5 <- mhglm_ml(y ~ x + (x | g1) + (x | g1:g2) + (x | g1:g2:g3),
                       data = g_data)
    model6 <- mhglm_ml(y ~ x + (x | g1:g2:g3) + (x | g1:g2) + (x | g1),
                       data = g_data)

    expect_equal(fixef(model1), fixef(model2))
    expect_equal(fixef(model1), fixef(model3))
    expect_equal(fixef(model1), fixef(model4))
    expect_equal(fixef(model1), fixef(model5))
    expect_equal(fixef(model1), fixef(model6))

    expect_equivalent(ranef(model1), ranef(model2))
    expect_equivalent(ranef(model1), ranef(model3))
    expect_equivalent(ranef(model1), ranef(model4))
    expect_equivalent(ranef(model1), ranef(model5))
    expect_equivalent(ranef(model1), ranef(model6))
})

test_that("test different colon permutations give the same random effects", {
    m1 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g2:g1) + (1|g3:g2:g1), data = g_data)
    m2 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g2:g1) + (1|g3:g1:g2), data = g_data)
    m3 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g2:g1) + (1|g1:g2:g3), data = g_data)
    m4 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g2:g1) + (1|g1:g3:g2), data = g_data)
    m5 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g2:g1) + (1|g2:g3:g1), data = g_data)
    m6 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g2:g1) + (1|g2:g1:g3), data = g_data)
    m7 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g1:g2) + (1|g3:g2:g1), data = g_data)
    m8 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g1:g2) + (1|g3:g1:g2), data = g_data)
    m9 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g1:g2) + (1|g1:g2:g3), data = g_data)
    m10 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g1:g2) + (1|g1:g3:g2), data = g_data)
    m11 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g1:g2) + (1|g2:g3:g1), data = g_data)
    m12 <- mhglm_ml(y ~ 1 + (1|g1) + (1|g1:g2) + (1|g2:g1:g3), data = g_data)

    expect_equivalent(ranef(m1), ranef(m2))
    expect_equivalent(ranef(m1), ranef(m3))
    expect_equivalent(ranef(m1), ranef(m4))
    expect_equivalent(ranef(m1), ranef(m5))
    expect_equivalent(ranef(m1), ranef(m6))
    expect_equivalent(ranef(m1), ranef(m7))
    expect_equivalent(ranef(m1), ranef(m8))
    expect_equivalent(ranef(m1), ranef(m9))
    expect_equivalent(ranef(m1), ranef(m9))
    expect_equivalent(ranef(m1), ranef(m10))
    expect_equivalent(ranef(m1), ranef(m11))
    expect_equivalent(ranef(m1), ranef(m12))
})

test_that("test order_bars function on properly specified formulas", {
    formulas <- c(
        y ~ 1 + (1|g1/g2/g3),

        y ~ 1 + (1|g1) + (1|g2:g1) + (1|g3:g2:g1),
        y ~ 1 + (1|g1) + (1|g3:g2:g1) + (1|g2:g1),
        y ~ 1 + (1|g2:g1) + (1|g1) + (1|g3:g2:g1),
        y ~ 1 + (1|g2:g1) + (1|g3:g2:g1) + (1|g1),
        y ~ 1 + (1|g3:g2:g1) + (1|g1) + (1|g2:g1),
        y ~ 1 + (1|g3:g2:g1) + (1|g2:g1) + (1|g1)
    )
    formula_bars <- Map(lme4::findbars, formulas)
    order_bars <- Map(order_bars_hierarchy, formula_bars)

    slash_result <- list(
        quote(1 | g1),
        quote(1 | g2:g1),
        quote(1 | g3:(g2:g1))
    )
    colon_result <- list(
        quote(1 | g1),
        quote(1 | g2:g1),
        quote(1 | g3:g2:g1)
    )

    expect_equal(order_bars[[1]], slash_result)

    expect_equal(order_bars[[2]], colon_result)
    expect_equal(order_bars[[3]], colon_result)
    expect_equal(order_bars[[4]], colon_result)
    expect_equal(order_bars[[5]], colon_result)
    expect_equal(order_bars[[6]], colon_result)
    expect_equal(order_bars[[7]], colon_result)
})

test_that("test order_bars function on properly misspecified formulas", {
    formulas <- c(
        y ~ 1 + (1|g1) + (1|g2) + (1|g3),
        y ~ 1 + (1|g1) + (1|g2:g1) + (1|g3:g1),
        y ~ 1 + (1|g4:g3:g2:g1) + (1|g2:g1) + (1|g1)
    )
    formula_bars <- Map(lme4::findbars, formulas)

    expect_error(order_bars_hierarchy(formula_bars[[1]]))
    expect_error(order_bars_hierarchy(formula_bars[[2]]))
    expect_error(order_bars_hierarchy(formula_bars[[3]]))
})

test_that("success using diagonal covariance", {
    model <- mhglm_ml(y ~ x + (x | g1/g2/g3), data = g_data,
                      control = mhglm_ml.control(diagcov = TRUE))
    varcor <- VarCorr(model)
    expect_equal(varcor[[1]][1, 2], 0)
    expect_equal(varcor[[2]][1, 2], 0)
    expect_equal(varcor[[3]][1, 2], 0)
})

test_that("succeeds on simdata for two levels, binomial", {
    model <- mhglm_ml(y ~ 1 + x + (1 + x | g1) + (1 + x | g2:g1), data = b_data)

    # fixef
    fixef0 <- c("(Intercept)" = 0.6351, "x" = 0.1447)
    fixef_m <- fixef(model)
    expect_equal_tol4(fixef_m, fixef0)

    # vcov
    vcov0 <- matrix(c(0.0016, -0.0004, -0.0004, 0.0021), nrow = 2L, ncol = 2L,
        dimnames = list(c("(Intercept)", "x"), c("(Intercept)", "x")))
    vcov_m <- vcov(model)
    expect_equal_tol4(vcov_m, vcov0)

    # VarCorr
    varcor0 <- list(
        g1 = matrix(c(0.0163, -0.0063, -0.0063, 0.0168), nrow = 2L, ncol = 2L,
            dimnames = list(c("(Intercept)", "x"), c("(Intercept)", "x"))),
        g2 = matrix(c(0.0250, -0.0004, -0.0004, 0.0301), nrow = 2L, ncol = 2L,
            dimnames = list(c("(Intercept)", "x"), c("(Intercept)", "x"))),
        sc = 0.41289)
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

    ranef0 <- structure(
        list(
            g1 = structure(
                list(
                    `(Intercept)` = c(
                        -0.0396, 0.1246, 0.0262, -0.1693, 0.089, 0.0842,
                        -0.0838, -0.1327, -0.082, 0.1788
                    ),
                    x = c(
                        -0.0041, -0.1768, -0.0069, -0.008, 0.115, -0.2117,
                        0.1244, 0.0548, 0.1201, -0.0036
                    )
                ),
                .Names = c("(Intercept)", "x"),
                row.names = as.character(1:10),
                class = "data.frame"
            ),
            `g2:g1` = structure(
                list(
                    `(Intercept)` = c(
                        0.2416, -0.1457, 0.0606, 2e-04, -0.1976, -0.16, 0.1418,
                        -0.1969, 0.0824, 0.0959, 0.0217, 0.0766, -0.2158,
                        0.3029, -0.1263, 0.0184, -0.0984, -0.1024, -0.2479,
                        0.1581, 0.1541, 0.1749, 0.0976, 0.074, -0.3265, -0.2264,
                        0.1108, -0.0665, 0.1274, -0.3439, -0.0941, 0.0546,
                        0.1199, -0.0293, -0.0336, -0.0081, 0.1346, -0.2952,
                        -0.1348, 0.0883, -0.1861, 0.0955, -0.1616, 0.1719,
                        0.0507, 0.1099, 0.1622, -0.1048, -0.052, -0.0106
                    ),
                    x = c(
                        -0.1203, 0.2411, -0.4445, 0.0249, 0.2719, 0.026, 0.0394,
                        -0.2696, 0.0645, 0.0464, -0.1775, 0.055, 0.1352,
                        -0.1063, 0.0541, -0.0333, 0.0837, -0.101, 0.0117,
                        0.0054, -0.0265, -0.1324, 0.0545, -0.0647, 0.1478,
                        -0.2373, 0.1438, 0.0038, -0.0434, -0.2918, 0.1443, 0.27,
                        0.1225, 0.2051, -0.3525, 0.0354, -0.0621, -0.2728,
                        0.1778, 0.0421, 0.0306, 0.1514, -0.0705, -0.024,
                        -0.0492, -0.039, -0.0854, -0.0083, 0.2272, -0.0875
                    )
                ),
                .Names = c("(Intercept)", "x"),
                row.names = c(
                    "1:1", "2:1", "3:1", "4:1", "5:1", "6:2", "7:2", "8:2",
                    "9:2", "10:2", "11:3", "12:3", "13:3", "14:3", "15:3",
                    "16:4", "17:4", "18:4", "19:4", "20:4", "21:5", "22:5",
                    "23:5", "24:5", "25:5", "26:6", "27:6", "28:6", "29:6",
                    "30:6", "31:7", "32:7", "33:7", "34:7", "35:7", "36:8",
                    "37:8", "38:8", "39:8", "40:8", "41:9", "42:9", "43:9",
                    "44:9", "45:9", "46:10", "47:10", "48:10", "49:10", "50:10"
                ),
                class = "data.frame")
            ),
        .Names = c("g1", "g2:g1"))

    ranef_m <- ranef(model)

    expect_equal(round(ranef_m$g1, 4), ranef0$g1)
    expect_equal(round(ranef_m$`g2:g1`, 4), ranef0$`g2:g1`)
})
