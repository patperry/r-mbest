
context("hglm.fast")

require(lme4)

test_that("succeeds on sleepstudy", {
    model <- hglm.fast(Reaction ~ Days, Subject, data=sleepstudy)
})
