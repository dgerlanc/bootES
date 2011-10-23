## Daniel Gerlanc and Kris Kirby (2010)
## Input and regression tests for the calcCohensD function

test.calcCohensD <- function() {
  ## Test the functioning of calcCohensD function

  ## Create 2-group data.frame.
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  twoGpsVec = c(g1, g2)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  freqVec   = rep(1, length(twoGpsVec))
  
  ## Regression test of calcCohensD
  truth = (mean(g1) - mean(g2)) / bootES:::pooledSD(twoGpsVec, grpLabels)
  d.res = bootES:::calcCohensD(twoGpsVec, freq=freqVec, grps=grpLabels)
    
  checkEquals(truth, d.res, tol=1e-4)

  ## Regression test of calcCohensD w/ cohens.d.sigma=TRUE

  ## Regression test of calcCohensD w/ glass.control=TRUE
  
}
