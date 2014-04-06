
context("hglm.fast")

require(lme4)

test_that("succeeds on sleepstudy", {
    model <- hglm.fast(Reaction ~ Days, Subject, data=sleepstudy)

    fixef0 <- c("(Intercept)" = 251.4, "Days" = 10.5)
    expect_that(round(fixef(model), 1), equals(fixef0))

    vcov0 <- matrix(c(44.0, -1.4, -1.4, 2.3), 2, 2)
    rownames(vcov0) <- colnames(vcov0) <- c("(Intercept)", "Days")
    expect_that(round(vcov(model), 1), equals(vcov0))

    varcor0 <- matrix(c(565.5, 11.1, 11.1, 32.7), 2, 2)
    rownames(varcor0) <- colnames(varcor0) <- c("(Intercept)", "Days")

    varcor <- VarCorr(model)[["Subject"]]
    expect_that(attr(varcor, "stddev"), equals(sqrt(diag(varcor))))

    expect_that(as.vector(attr(varcor, "correlation")),
                equals(as.vector(t(varcor / attr(varcor, "stddev"))
                                 / attr(varcor, "stddev"))))

    attr(varcor, "stddev") <- NULL
    attr(varcor, "correlation") <- NULL
    expect_that(round(varcor, 1), equals(varcor0))
})
