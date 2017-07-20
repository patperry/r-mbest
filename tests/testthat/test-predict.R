context("mhglm-parallel")

test_that("lme4::sleepstudy predictions regression test", {
    m_seq <- mhglm(Reaction ~ Days + (Days | Subject), data=lme4::sleepstudy,
                 control=list(parallel=FALSE))
    #basicConfig("INFO")
    m_par <- mhglm(Reaction ~ Days + (Days | Subject), data=lme4::sleepstudy,
                 control=list(parallel=TRUE))
    
    sequential_predictions <- predict(m_seq, lme4::sleepstudy, se.fit=TRUE)
    parallel_predictions <- predict(m_par, lme4::sleepstudy, se.fit=TRUE)
    
    expect_equal(sequential_predictions[['fit']], 
                 parallel_predictions[['fit']])
    expect_equal(sequential_predictions[['se.fit']], 
                 parallel_predictions[['se.fit']])
    expect_equal(sequential_predictions[['residual.scale']], 
                 parallel_predictions[['residual.scale']])
})


test_that("Simulated data predictions regression test", {
    # example taken from http://www.r-bloggers.com/random-regression-coefficients-using-lme4/
    set.seed(5432)
    J <- 20
    N <- 3000
    train.df <- data.frame( unit = sort(rep(c(1:N),J)), 
                           J = rep(c(1:J),N) , x = rnorm(n = J*N))
    beta <- 3 + .2*rnorm(N)
    train.df$beta <- beta[train.df$unit]
    train.df$y <- 1 + train.df$x * train.df$beta + .75*rnorm(n = J*N)
    
    test.df <- train.df
    test.df$y <- 1 + test.df$x * test.df$beta + .75*rnorm(n = J*N)
  
    suppressWarnings({
        m_seq <- mhglm(y ~ 1 + x + (1+x|unit), data = train.df,
                       control=list(parallel=FALSE))
    })

    #basicConfig("INFO")
    suppressWarnings({
        m_par <- mhglm(y ~ 1 + x + (1 + x | unit), data=train.df,
                       control=list(parallel=TRUE))
    })
    
    sequential_predictions <- predict(m_seq, test.df, se.fit=TRUE)
    parallel_predictions <- predict(m_par, test.df, se.fit=TRUE)
    
    expect_equal(sequential_predictions[['fit']], 
                 parallel_predictions[['fit']])
    expect_equal(sequential_predictions[['se.fit']], 
                 parallel_predictions[['se.fit']])
    expect_equal(sequential_predictions[['residual.scale']], 
                 parallel_predictions[['residual.scale']])
})

