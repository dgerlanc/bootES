## Daniel Gerlanc and Kris Kirby (2010)
## Input and regression tests for the bootES function

test.bootES.input <- function() {
  ## Test the functioning of the bootES interface w/ invalid inputs.

  g1 = c(11, 12, 13, 14, 15)
  g2 = c(26, 27, 28, 29)
  g3 = c(17, 18, 19, 20, 21, 22, 23)
  
  grpLabels = rep(c("A", "B", "C"), times=c(length(g1), length(g2), length(g3)))
  threeGps  = data.frame(grpLabels, scores=c(g1, g2, g3))
  threeGpsVec = c(g1, g2, g3)
  lambdas   = c(A=-1, B=2, C=-1)
    
  ## Pass a non-data.frame object as 'dat'
  res = try(bootES("foo"), silent=TRUE)
  checkTrue(grepl("'dat' must be a data.frame or numeric vector.", res))
  
  ## Pass a data.frame to 'dat' with no records
  checkException(bootES(data.frame()), silent=TRUE)
  
  ## Pass an 'R' of length greater than 1
  checkException(bootES(data.frame(scores=1), R=c(2, 1)), silent=TRUE)
  
  ## Pass an 'R' with a fractional portion
  res <- try(bootES(data.frame(scores=1), R=0.5), silent=TRUE)
  checkTrue(grepl("integer of length 1", res))
  
  ## Use a 'data.col' not in 'dat'
  checkException(bootES(threeGps, R=100, data.col="foo"), silent=TRUE)

  ## Use a 'glass.control' value that is not a valid group
  res <- try(bootES(data.frame(scores=1), glass.control="foo"), silent=TRUE)
  checkTrue(grepl("'glass.control' missing", res))
}

test.bootES.univariate <- function() {
  ## Test the functioning of the univariate bootstrap functions through the
  ## bootES interface.

  g1 = c(11, 12, 13, 14, 15)
  g2 = c(26, 27, 28, 29)
  g3 = c(17, 18, 19, 20, 21, 22, 23)
  
  grpLabels = rep(c("A", "B", "C"), times=c(length(g1), length(g2), length(g3)))
  threeGps  = data.frame(grpLabels, scores=c(g1, g2, g3))
  threeGpsVec = c(g1, g2, g3)
  lambdas   = c(A=-1, B=2, C=-1)
  
  ## Test: 'meanBoot' through 'bootES'
  set.seed(1)
  truth    = mean(threeGps$scores)
  mean.res = bootES(threeGps, R=1000, data.col="scores", effect.type="unstandardized")
  mean.res.vec = bootES(threeGpsVec, effect.type="unstandardized")
  checkEquals(truth, mean.res$t0)
  checkEquals(truth, mean.res.vec$t0)

  ## Test: 'rMeanBoot' through 'bootES'
  set.seed(1)
  truth     = bootES:::rMean(threeGps$scores)
  rMean.res = bootES(threeGps, R=1000, data.col="scores", effect.type="r")
  checkEquals(truth, rMean.res$t0)

  ## Test: 'dMeanBoot' through 'bootES'
  set.seed(1)
  truth     = bootES:::dMean(threeGps$scores)
  rMean.res = bootES(threeGps, R=1000, data.col="scores", effect.type="cohens.d")
  checkEquals(truth, rMean.res$t0)

  ## Test: 'dMeanBoot' and Cohen's Sigma d through 'bootES'
  set.seed(1)
  truth     = bootES:::dSigmaMeanBoot(threeGps$scores, 1:length(threeGps$scores))
  rMean.res = bootES(threeGps, R=1000, data.col="scores", effect.type="cohens.d.sigma")
  checkEquals(truth, rMean.res$t0)
  
  ## Test: 'dMeanBoot' and Hedge's g through 'bootES'
  set.seed(1)
  truth     = bootES:::dMean(threeGps$scores)
  rMean.res = bootES(threeGps, R=1000, data.col="scores", effect.type="hedges.g")
  checkEquals(truth, rMean.res$t0)
    
}

test.bootES.verbosity <- function() {
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  twoGpsA   = data.frame(x=c(g1, g2), team=grpLabels)
  lambdas   = c(A=1, B=-1)
  
  ## Integration test of stat='contrast' and effect.type='unstandardized'
  set.seed(1)
  truth     = mean(g1) - mean(g2)
  unstdDiff.res = bootES(twoGpsA, R=1000, data.col="x", grp.col="team",
    stat="contrast", effect.type="unstandardized", contrasts=lambdas)
  
}

test.bootES.multivariate <- function() {
  ## Test the functioning of the multivariate bootstrap functions through the
  ## bootES interface.

  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  twoGpsA   = data.frame(x=c(g1, g2), team=grpLabels)
  twoGpsErr = data.frame(x=c(g1, g2), team=rep("A", length(c(g1, g2))))
  lambdas   = c(A=1, B=-1)
  
  ## Integration test of stat='contrast' and effect.type='unstandardized'
  set.seed(1)
  truth     = mean(g1) - mean(g2)
  unstdDiff.res = bootES(twoGpsA, R=1000, data.col="x", grp.col="team",
    stat="contrast", effect.type="unstandardized", contrasts=lambdas)
  checkEquals(truth, unstdDiff.res$t0)  
  
  ## Integration test of stat='contrast' and effect.type='unstandardized' where
  ## there is only one group. This should cause an error.
  unstdDiff.err = try(bootES(twoGpsErr, R=1000, data.col="x", grp.col="team",
    stat="contrast", effect.type="unstandardized"), silent=TRUE)
  
  ## Integration test of stat='contrast' and effect.type='r'
  ## ToDo!
  path = system.file("gender.csv", package="bootES")
  gender = read.csv(path, strip.white=TRUE, header=TRUE)
  
  ## Integration test of stat='cor'
  set.seed(1)
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  twoGps    = data.frame(g1=g1)
  twoGps$g2 = rep(g2, length.out=nrow(twoGps))
  truth     = with(twoGps, cor(g1, g2))
  cor.res   = bootES(twoGps, R=10, stat="cor")
  checkEquals(truth, cor.res$t0)

  ## Integration test of stat='cor.diff'
  set.seed(1)
  a1 = c(1:5, -(1:5))
  a2 = c(10:14, 10:14)
  twoGps       = data.frame(a1, a2, group=rep(c(1, 2), each=5))
  truth        = cor(a1[1:5], a2[1:5]) - cor(a1[6:10], a2[6:10])
  cor.diff.res = bootES(twoGps, R=10, stat="cor.diff", grp.col="group")
  checkEquals(truth, cor.diff.res$t0)
}

test.bootES.mean <- function() {
  ## Integration tests of stat='mean' w/ other effect.types
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  twoGpsA   = data.frame(x=c(g1, g2), team=grpLabels)
  twoGpsErr = data.frame(x=c(g1, g2), team=rep("A", length(c(g1, g2))))
  lambdas   = c(A=1, B=-1)
  
  effect.types = eval(formals(bootES)[["effect.type"]])
  for (type in effect.types) {
    next
  }
}

test.bootES.verbose <- function() {
  ## Test the verbosity of bootES functions
  g1 = c(11, 12, 13, 14, 15)
  g2 = c(26, 27, 28, 29)
  g3 = c(17, 18, 19, 20, 21, 22, 23)
  
  grpLabels = rep(c("A", "B", "C"), times=c(length(g1), length(g2), length(g3)))
  threeGps  = data.frame(grpLabels, scores=c(g1, g2, g3))
  threeGpsVec = c(g1, g2, g3)
  lambdas   = c(A=-1, B=2, C=-1)

  ## It shouldn't print anything to the screen.
  test = capture.output(assign("mean.res", bootES(threeGps,
    data.col="scores", effect.type="unstandardized", verbose=0)))

  checkTrue(length(test) == 0)
  
  ## It should print the statistic and the confidence interval to the screen
  test = capture.output(assign("mean.res", bootES(threeGps,
    data.col="scores", effect.type="unstandardized", verbose=1)))
  
  checkTrue(grepl(" +Stat +CI \\(Low\\) +CI \\(High\\)", test[2], perl=TRUE))
  checkTrue(grepl("\\[1,\\] [\\d.]+ +[\\d.]+ +[\\d.]+", test[3], perl=TRUE))
}

test.bootES.automagic <- function() {
    return
}

test.meanUnweightedBoot <- function() {
  options(warn=2)
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  twoGpsA   = data.frame(x=c(g1, g2), team=grpLabels)

  bootES:::meanUnweightedBoot(twoGpsA$x, freq=rep(1, nrow(twoGpsA)),
                              grps=grpLabels)
}

test.determineStat <- function() {
  ## Test the function that determines which statistic to use
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  twoGpsA   = data.frame(x=c(g1, g2), team=grpLabels)
  twoGpsErr = data.frame(x=c(g1, g2), team=rep("A", length(c(g1, g2))))
  lambdas   = c(A=1, B=-1)

  ## Stat should be mean
  test = data.frame(score=c(g1, g2))
  checkEquals(bootES:::determineStat(test), 'mean')
  
  ## Stat should be slope
  test = data.frame(x=c(g1, g2), y=c(-g1, -g2))
  checkEquals(bootES:::determineStat(test, effect.type='slope'), 'slope')
  
  ## Stat should be cor
  test = data.frame(x=c(g1, g2), y=c(-g1, -g2))
  checkEquals(bootES:::determineStat(test), 'cor')

  ## Stat should be contrast
  test.dat = data.frame(score=c(g1, g2), group=grpLabels)
  test = bootES:::determineStat(test, grps=grpLabels, contrasts=lambdas)
  checkEquals(test, 'contrast')
  
  ## Stat should be cor.diff
  test.dat = data.frame(score=c(g1, g2), group=grpLabels)
  test = bootES:::determineStat(test, grps=grpLabels)
  checkEquals(test, 'cor.diff')
}
