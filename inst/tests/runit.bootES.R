## Daniel Gerlanc and Kris Kirby (2010)
## Input and regression tests for the bootES function

test.AAA <- function() {
  gender <<- read.csv(system.file("gender.csv", package="bootES"),
                      strip.white=TRUE, header=TRUE)
}

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

  ## Pass an invalid 'group' to 'contrasts'
  res = try(bootES(gender, data.col="Meas3",
    group.col="Condition", contrasts = c(Fake = -50, C = 50),
    scale.weights=TRUE), silent=TRUE)
  checkTrue(grepl("'Fake' is/are not valid groups.", res[1]))
  
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
  mean.res = bootES(threeGps, R=1000, data.col="scores",
    effect.type="unstandardized")
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
  rMean.res = bootES(threeGps, R=1000, data.col="scores",
    effect.type="cohens.d")
  checkEquals(truth, rMean.res$t0)

  ## Test: 'dMeanBoot' and Cohen's Sigma d through 'bootES'
  set.seed(1)
  truth     = bootES:::dSigmaMeanBoot(threeGps$scores,
    1:length(threeGps$scores))
  rMean.res = bootES(threeGps, R=1000, data.col="scores",
    effect.type="cohens.d.sigma")
  checkEquals(truth, rMean.res$t0)
  
  ## Test: 'dMeanBoot' and Hedge's g through 'bootES'
  set.seed(1)
  truth     = bootES:::dMean(threeGps$scores)
  rMean.res = bootES(threeGps, R=1000, data.col="scores",
    effect.type="hedges.g")
  checkEquals(truth, rMean.res$t0)
    
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
  unstdDiff.res = bootES(twoGpsA, R=1000, data.col="x", group.col="team",
    effect.type="unstandardized", contrasts=lambdas)
  checkEquals(truth, unstdDiff.res$t0)  
  
  ## Integration test of stat='contrast' and effect.type='unstandardized' where
  ## there is only one group. This should cause an error.
  unstdDiff.err = try(bootES(twoGpsErr, R=1000, data.col="x", group.col="team",
    effect.type="unstandardized"), silent=TRUE)    
  
  ## Integration test of stat='cor'
  set.seed(1)
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  twoGps    = data.frame(g1=g1)
  twoGps$g2 = rep(g2, length.out=nrow(twoGps))
  truth     = with(twoGps, cor(g1, g2))
  cor.res   = suppressWarnings(bootES(twoGps, R=10, effect.type="r"))
  checkEquals(truth, cor.res$t0)
  
}

test.bootES.contrast <- function() {

  ## Assert: Calculated value matches known value for an unstandardized
  ## contrast
  set.seed(1)
  truth = -522.43
  test  = bootES(gender, data.col="Meas3", group.col="Condition",
    contrasts = c(A = -40, B = -10, C = 50), scale.weights=FALSE)
  checkEquals(truth, test$t0, tol=1e-2)

  ## Assert: Calculated value matches known value for an unstandardized
  ## contrast with weights scaled
  set.seed(1)
  truth.contrast.scaled = -10.4486
  test  = bootES(gender, data.col="Meas3", group.col="Condition",
    contrasts = c(A = -40, B = -10, C = 50), scale.weights=TRUE)
  checkEquals(truth.contrast.scaled, test$t0, tol=1e-4)

  ## Assert: Calculated value matches known value for an unstandardized
  ## contrast with weights scaled and a group left out
  set.seed(1)
  truth.contrast.omit = -3.0535
  test  = bootES(gender, data.col="Meas3", group.col="Condition",
    contrasts = c(A = -1, C = 1))  
  checkEquals(truth.contrast.omit, test$t0, tol=1e-4)

  ## Assert: Default weights of -1 and 1 are used when not passed in
  test.dflt  = bootES(gender, data.col="Meas3", group.col="Condition",
    contrasts = c('A', 'C'))
  checkEquals(truth.contrast.omit, test.dflt$t0, tol=1e-4)
}

test.bootES.cor.diff <- function() {
  ## Integration test of stat='cor.diff'
  ## Note: Tests variable ordering of numberic and group columns.
  set.seed(1)
  a1 = c(1:5, -(1:5))
  a2 = c(10:14, 10:14)
  twoGps       = data.frame(group=rep(c(1, 2), each=5), a1, a2)
  truth        = cor(a1[1:5], a2[1:5]) - cor(a1[6:10], a2[6:10])
  cor.diff.res = bootES(twoGps, R=10, group.col="group", effect.type="r")
  checkEquals(truth, cor.diff.res$t0)
}

test.bootES.slope <- function() {
  ## Regression test for when effect.type='slope'
  set.seed(1)
  truth <- -0.1244
  test  <- bootES(gender, data.col="Meas3", group.col="Condition",
                  effect.type="slope",
                  slope.levels=c(A=30, B=60, C=120))
  checkEquals(truth, test$t0, tol=1e-2)
}

test.bootES.output <- function() {
  
  ## Test the verbosity of bootES functions
  g1 = c(11, 12, 13, 14, 15)
  g2 = c(26, 27, 28, 29)
  g3 = c(17, 18, 19, 20, 21, 22, 23)
  
  grpLabels = rep(c("A", "B", "C"), times=c(length(g1), length(g2), length(g3)))
  threeGps  = data.frame(grpLabels, scores=c(g1, g2, g3))
  threeGpsVec = c(g1, g2, g3)
  lambdas   = c(A=-1, B=0, C=1)

  ## It shouldn't print anything to the screen.
  test = capture.output(assign("mean.res", bootES(threeGps,
    data.col="scores", effect.type="unstandardized")))

  checkTrue(length(test) == 0)
  
  ## Validate output
  test = capture.output(print(bootES(threeGps,
    data.col="scores", effect.type="unstandardized")))
  
  str1 = " +Stat +CI \\(Low\\) +CI \\(High\\) +bias +std\\. error"
  checkTrue(grepl(str1, test[2], perl=TRUE))
  str2 = "\\[1,\\] [\\d.]+ +[\\d.]+ +[\\d.]+ +[-\\d.]+ +[\\d.]+"
  checkTrue(grepl(str2, test[3], perl=TRUE))

  ## It should print the statistic and the confidence interval to the screen
  test = capture.output(print(bootES(gender, data.col="Meas3",
    group.col="Condition", contrasts = c(A = -40, B = -10, C = 50),
    scale.weights=TRUE)))
  
  checkTrue(grepl("User-specified lambdas: \\(-?\\d+, -?\\d+, -?\\d+\\)", 
                  test[1], perl=TRUE))
  regexp = paste("Scaled lambdas: \\(-?[0-9.]+, -?[0-9.]+, -?[0-9.]+\\)")
  checkTrue(grepl(regexp, test[2], perl=TRUE))
}

test.bootES.citype <- function() {
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
  mean.res = bootES(threeGps, R=1000, data.col="scores",
    effect.type="unstandardized")
  mean.res.vec = bootES(threeGpsVec, effect.type="unstandardized")
  checkEquals(truth, mean.res$t0)
  checkEquals(truth, mean.res.vec$t0)

  ci.types = eval(formals(bootES)$ci.type)
  ci.types = ci.types[!ci.types %in% "stud"]
  for (ci.type in ci.types) {
    set.seed(1)
    if (ci.type == "stud") {
      . <- bootES(threeGps, R=1000, data.col="scores",
                  effect.type="unstandardized", ci.type=ci.type,
                  var.t0=1, var.t=1)
    } else {
      . <- bootES(threeGps, R=1000, data.col="scores",
                  effect.type="unstandardized", ci.type=ci.type)
    }
  }
}
