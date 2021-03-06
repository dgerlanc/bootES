## Daniel Gerlanc and Kris Kirby (2010-2021)

test_that("determineStat calculates correct statistic to use", {
  ## Test the function that determines which statistic to use
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  g3        = c(17, 18, 19, 20, 21, 22, 23)
  grpLabels      = rep(c("A", "B"), times=c(length(g1), length(g2)))
  threeGrpLabels = rep(c("A", "B", "C"),
    times=c(length(g1), length(g2), length(g3)))
  lambdas   = c(A=1, B=-1)

  ## 1.1.r One Group, One Measure on One Dependent Variable
  ## Assert: stat is mean
  test.input  = data.frame(score=c(g1, g2))
  test.result <- bootES:::determineStat(test.input,
                                        effect.type="unstandardized")
  expect_equal(test.result, 'mean')

  ## Assert that Three Groups, One Measure, Unstandardized yields mean
  ## Assert: stat is mean
  test.input = data.frame(score=c(g1, g2, g3), group=threeGrpLabels)
  test.result <- bootES:::determineStat(test.input$score,
                                        grps=test.input$group,
                                        effect.type="unstandardized")
  expect_equal(test.result, 'mean')

  ## 1.2.r One Group, Two Measures on Two Different Dependent Variables
  ## Assert: stat is slope
  test.input  <- data.frame(x=c(g1, g2), y=c(-g1, -g2))
  test.result <- bootES:::determineStat(test.input,
                                        data.col=test.input$x,
                                        grps=test.input$y,
                                        effect.type='slope')
  expect_equal(test.result, 'slope')
      
  ## Assert: stat is contrast
  test.input = data.frame(score=c(g1, g2), group=grpLabels)
  test.result <- bootES:::determineStat(test.input,
                                        grps=grpLabels,
                                        contrast=lambdas)
  expect_equal(test.result, 'contrast')
    
})

test_that("correct statistic used for correlation", {
  g1 = c(11, 12, 13, 14, 15)
  g2 = c(26, 27, 28, 29)
  
  ## Assert: stat is cor when effect.type 'r' selected
  test.input  <- data.frame(x=c(g1, g2), y=c(-g1, -g2))
  test.result <- bootES:::determineStat(test.input,
                                        effect.type="r")
  expect_equal(test.result, 'cor')

  ## Assert: stat is cor when effect.type 'unstandardized' is selected
  test.input  <- data.frame(x=c(g1, g2), y=c(-g1, -g2))
  test.result <- bootES:::determineStat(test.input,
                                        effect.type="unstandardized")
  expect_equal(test.result, 'cor')
})

test_that("correct statistic used for cor.diff", {
  g1 = c(11, 12, 13, 14, 15)
  g2 = c(26, 27, 28, 29)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  
  ## Assert: stat is cor.diff
  test.input <- data.frame(x=c(g1, g2), y=c(-g1, -g2), group=grpLabels)
  test.result <- bootES:::determineStat(test.input,
                                        grps=test.input$group,
                                        effect.type="r")
  expect_equal(test.result, 'cor.diff')
})
