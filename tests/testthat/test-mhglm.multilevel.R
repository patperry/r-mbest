context("mhglm")
library("lme4")

test_that("succeeds on sleepstudy with mhglm_ml", {
    model <- mhglm_ml(Reaction ~ Days + (Days | Subject), data=sleepstudy)

    # fixef
    fixef0 <- c("(Intercept)" = 251.4, "Days" = 10.5)
    expect_that(round(fixef(model), 1), equals(fixef0))

    # vcov
    vcov0 <- matrix(c(44.0, -1.4, -1.4, 2.3), 2, 2)
    rownames(vcov0) <- colnames(vcov0) <- c("(Intercept)", "Days")
    expect_that(round(vcov(model), 1), equals(vcov0))

    # VarCorr
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

    # ranef
    ranef0 <- matrix(c( 2.8, -40.0, -38.4, 22.8, 21.5,  8.8,  16.4, -7.0,  -1.0,
                       34.7, -24.6, -12.3,  4.3, 20.6,  3.3, -24.7,  0.7,  12.1,
                        9.1,  -8.6,  -5.5, -4.7, -2.9, -0.2,  -0.2,  1.0, -10.6,
                        8.6,   1.1,   6.5, -3.0,  3.6,  0.9,   4.7, -1.0,   1.3),
                     18, 2)
    rownames(ranef0) <- as.character(c(308, 309, 310, 330, 331, 332, 333, 334,
                                       335, 337, 349, 350, 351, 352, 369, 370,
                                       371, 372))
    colnames(ranef0) <- c("(Intercept)", "Days")
#    expect_that(round(as.matrix(ranef(model)[["Subject"]]), 1), equals(ranef0)) #TODO: need to add column name
    expect_that(round(as.matrix(ranef(model)[[1]]), 1), equals(ranef0))
})



test_that("succeeds on sleepstudy with two hiererchies", {
    set.seed(0)
    fakegroup <- paste0(sleepstudy[,'Subject'],"_",sample(c(1:2),180,replace = TRUE))
    sleepstudy2 <- data.frame(sleepstudy, fakegroup)
    model <- suppressWarnings(mhglm_ml(Reaction ~ Days + (Days | Subject ) + (Days |  fakegroup), data = sleepstudy2,
                      control = list(standardize = FALSE, diagcov = FALSE)))


    # fixef
    fixef0 <- c("(Intercept)"= 251.4, "Days" = 10.3)
    expect_that(round(fixef(model),1), equals(fixef0))

    # vcov
    vcov0 <- matrix(c(67.0, 1.3, 1.3, 2.0), 2, 2)
    rownames(vcov0) <- colnames(vcov0) <- c("(Intercept)","Days")
    expect_that(round(vcov(model), 1), equals(vcov0))

    # VarCorr
    varcor0_1 <- matrix(c(653.4, 7.9, 7.9, 32.7), 2, 2)
    rownames(varcor0_1) <- colnames(varcor0_1) <- c("(Intercept)","Days")
    varcor0_2 <- matrix(c(32.8, -5.5, -5.5, 1.0), 2, 2)
    rownames(varcor0_2) <- colnames(varcor0_2) <- c("(Intercept)","Days")
    varcor1 <- VarCorr(model)[[1]]
    varcor2 <- VarCorr(model)[[2]]
    attr(varcor1,'stddev')<-NULL; attr(varcor1,'correlation')<-NULL
    attr(varcor2,'stddev')<-NULL; attr(varcor2,'correlation')<-NULL

    expect_that(round(varcor1,1), equals(varcor0_1))
    expect_that(round(varcor2,1), equals(varcor0_2))


    # ranef
    ranef.est  <- ranef(model)
    paste0(round(c(as.matrix(ranef.est[[1]])),1),collapse = ', ')
    paste0(round(c(as.matrix(ranef.est[[2]])),1),collapse = ', ')
    ranef0_1 <- matrix(c(9.8, -45.4, -45.1, 22.4, 23, 1.1, 21.4, -8, 5.1, 37, 
                         -31, -13.4, 5.8, 22.5, 3.5, -31.3, 0.9, 14.8, 7, -8.1, 
                         -4.6, -4.6, -3, 0.8, -0.8, 1.3, -12.1, 8.6, 2.3, 6.8, 
                         -3.1, 3.5, 1, 6, -0.9, 1.1),ncol = 2)
    rownames(ranef0_1) <- as.character(c(308, 309, 310, 330, 331, 332, 333, 334,
                                       335, 337, 349, 350, 351, 352, 369, 370,
                                       371, 372))
    colnames(ranef0_1) <- c("(Intercept)","Days")

    ranef0_2 <- matrix(c(-0.5, -0.3, 0.2, -1.1, 0.7, -1.6, 4.9, -1.5, 1.2, 0.4, 
                         0, 0.5, -1.1, 2.4, 0, -0.6, 8.9, -3.2, 0.1, 0.2, 
                         -1.3, -0.5, -1.7, -0.1, -1.4, 2.1, 1.9, -1.4, -0.8, 0.7, 
                         -2.4, -0.1, 1, -0.8, -0.7, 1.1, 0.1, 0.1, 0, 0.2, 
                         -0.1, 0.2, -0.9, 0.3, -0.2, 0, 0, -0.1, 0.2, -0.4, 
                         0, 0.1, -1.5, 0.5, 0, 0, 0.2, 0.1, 0.3, 0, 
                         0.3, -0.4, -0.3, 0.2, 0.1, -0.1, 0.4, 0.1, -0.2, 0.2, 0.1, -0.1),ncol = 2)
    rownames(ranef0_2) <- c("308_1", "308_2", "309_1", "309_2", "310_1", "310_2", "330_1", "330_2", "331_1", 
                          "331_2", "332_1", "332_2", "333_1", "333_2", "334_1", "334_2", "335_1", "335_2", 
                          "337_1", "337_2", "349_1", "349_2", "350_1", "350_2", "351_1", "351_2", "352_1",
                          "352_2", "369_1", "369_2", "370_1", "370_2", "371_1", "371_2", "372_1", "372_2")
    colnames(ranef0_2) <- c("(Intercept)","Days")

    expect_that(round(as.matrix(ranef.est[[1]]),1), equals(ranef0_1))
    expect_that(round(as.matrix(ranef.est[[2]]),1), equals(ranef0_2))

})

test_that("succeeds on parsing three-level models, using different formulas",{
    set.seed(0)
    fakegroup <- sample(c(1:2),180,replace = TRUE)
    fakegroup2 <- paste0(sleepstudy[,'Subject'],"_",fakegroup)
    sleepstudy2 <- data.frame(sleepstudy, fakegroup, fakegroup2)
    control <- list(standardize = FALSE, diagcov = FALSE)
    model1 <- suppressWarnings(mhglm_ml(Reaction ~ Days + (Days | Subject) + (Days | fakegroup2 ), data = sleepstudy2, control = control))
    model2 <- suppressWarnings(mhglm_ml(Reaction ~ Days + (Days | Subject/fakegroup ), data = sleepstudy2, control = control))
    model3 <- suppressWarnings(mhglm_ml(Reaction ~ Days + (Days | Subject/fakegroup2 ), data = sleepstudy2,control = control))
    model4 <- suppressWarnings(mhglm_ml(Reaction ~ Days + (Days | Subject) + (Days|Subject:fakegroup ), data = sleepstudy2,control = control))
    model5 <- suppressWarnings(mhglm_ml(Reaction ~ Days + (Days | Subject) + (Days|Subject:fakegroup2 ), data = sleepstudy2,control = control))


    # fixef
    fixef0 <- c("(Intercept)"= 251.4, "Days" = 10.3)
    expect_that(round(fixef(model1),1), equals(fixef0))
    expect_that(round(fixef(model2),1), equals(fixef0))
    expect_that(round(fixef(model3),1), equals(fixef0))
    expect_that(round(fixef(model4),1), equals(fixef0))
    expect_that(round(fixef(model5),1), equals(fixef0))

    # ranef
    ranef.est1  <- ranef(model1)
    ranef.est2  <- ranef(model2)
    ranef.est3  <- ranef(model3)
    ranef.est4  <- ranef(model4)
    ranef.est5  <- ranef(model5)
    ranef0_1 <- matrix(c(9.8, -45.4, -45.1, 22.4, 23, 1.1, 21.4, -8, 5.1, 37, 
                         -31, -13.4, 5.8, 22.5, 3.5, -31.3, 0.9, 14.8, 7, -8.1, 
                         -4.6, -4.6, -3, 0.8, -0.8, 1.3, -12.1, 8.6, 2.3, 6.8, 
                         -3.1, 3.5, 1, 6, -0.9, 1.1),ncol = 2)
    rownames(ranef0_1) <- as.character(c(308, 309, 310, 330, 331, 332, 333, 334,
                                       335, 337, 349, 350, 351, 352, 369, 370,
                                       371, 372))
    colnames(ranef0_1) <- c("(Intercept)","Days")

    ranef0_2 <- matrix(c(-0.5, -0.3, 0.2, -1.1, 0.7, -1.6, 4.9, -1.5, 1.2, 0.4, 
                         0, 0.5, -1.1, 2.4, 0, -0.6, 8.9, -3.2, 0.1, 0.2, 
                         -1.3, -0.5, -1.7, -0.1, -1.4, 2.1, 1.9, -1.4, -0.8, 0.7, 
                         -2.4, -0.1, 1, -0.8, -0.7, 1.1, 0.1, 0.1, 0, 0.2, 
                         -0.1, 0.2, -0.9, 0.3, -0.2, 0, 0, -0.1, 0.2, -0.4, 
                         0, 0.1, -1.5, 0.5, 0, 0, 0.2, 0.1, 0.3, 0, 
                         0.3, -0.4, -0.3, 0.2, 0.1, -0.1, 0.4, 0.1, -0.2, 0.2, 0.1, -0.1),ncol = 2)
    colnames(ranef0_2) <- c("(Intercept)","Days")

    expect_that(round(as.matrix(ranef.est1[[1]]),1), equals(ranef0_1))
    expect_that(round(as.matrix(ranef.est2[[1]]),1), equals(ranef0_1))
    expect_that(round(as.matrix(ranef.est3[[1]]),1), equals(ranef0_1))
    expect_that(round(as.matrix(ranef.est4[[1]]),1), equals(ranef0_1))
    expect_that(round(as.matrix(ranef.est5[[1]]),1), equals(ranef0_1))

    tmp <- round(as.matrix(ranef.est1[[2]]),1); rownames(tmp) <- NULL; expect_that(tmp, equals(ranef0_2))
    tmp <- round(as.matrix(ranef.est3[[2]]),1); rownames(tmp) <- NULL; expect_that(tmp, equals(ranef0_2))
    tmp <- round(as.matrix(ranef.est4[[2]]),1); rownames(tmp) <- NULL; expect_that(tmp, equals(ranef0_2))
    tmp <- round(as.matrix(ranef.est5[[2]]),1); rownames(tmp) <- NULL; expect_that(tmp, equals(ranef0_2))

    # reverse the ranef group name from fakegroup:Study to Study:fakegroup, then sort by rownames
    tmp <- round(as.matrix(ranef.est2[[2]]),1); 
    rownames(tmp) <- apply(array(rownames(tmp)) , 1, function(x){
                               y = strsplit(x,":")[[1]];
                               paste0(y[2] ,":",y[1]) })
    tmp <- tmp[order(rownames(tmp)),]; rownames(tmp) <- NULL; expect_that(tmp, equals(ranef0_2)) 
})

