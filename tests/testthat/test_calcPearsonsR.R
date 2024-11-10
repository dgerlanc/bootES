## Daniel Gerlanc and Kris Kirby (2010-2012)
## Input and regression tests for the calcPearsonsR function

test_that("calcPearsonsR produces known result", {

  ## Create 2-group data.frame.
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  g3        = c(17, 18, 19, 20, 21, 22, 23)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  twoGpsA   = data.frame(x=c(g1, g2), team=grpLabels)
  twoGpsErr = data.frame(x=c(g1, g2), team=rep("A", length(c(g1, g2))))
  lambdas   = c(A=1, B=-1)

  ## Create 3-group data.frame.
  grpLabels3  = rep(c("A", "B", "C"), times=c(length(g1), length(g2), length(g3)))
  threeGps    = data.frame(grpLabels3, scores=c(g1, g2, g3))
  threeGpsVec = c(g1, g2, g3)
  lambdas3    = c(A=-1, B=2, C=-1)
  
  ## Regression test of calcPearsonsR for a 3-group contrast
  set.seed(1)
  truth = 0.9452
  r.res = bootES:::calcPearsonsR(threeGpsVec, freq=rep(1, length(threeGpsVec)), grps=grpLabels3, lambdas=lambdas3)
  expect_equal(truth, r.res, tolerance=1e-4)
})

test_that("produces same result with non-lexicographic contrast sorting", {
  dat = data.frame(cond=rep(c("A", "B", "C"), each=10), score=1:30)

  truth = 0.926
  res_1 = bootES(
    dat, data.col="score", group.col="cond", contrast=c(A=-1, B=0.5, C=0.5),
    effect.type="r")$t0
  expect_equal(truth, res_1, tolerance=1e-3)


  res_2 = bootES(
    dat, data.col="score", group.col="cond", contrast=c(B=0.5, C=0.5, A=-1),
    effect.type="r")$t0
  expect_equal(truth, res_2, tolerance=1e-3)


  res_3 = bootES(
    dat, data.col="score", group.col="cond", contrast=c(C=0.5, A=-1, B=0.5),
    effect.type="r")$t0
  expect_equal(truth, res_3, tolerance=1e-3)
})
