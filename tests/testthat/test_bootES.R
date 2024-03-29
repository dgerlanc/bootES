## Daniel Gerlanc and Kris Kirby (2010-2012)
## Input and regression tests for the bootES function

get_gender <- function() {
  gender <- read.csv(system.file("example.csv", package="bootES"),
                      strip.white=TRUE, header=TRUE)
  gender$GenderByCond <- paste(gender$Gender, gender$Condition, sep = "-")
  gender
}

test_that("bootES handles invalid 'R' values", {
  r_check_message = paste0(
    "R must be of 'type' integer or 'type' real with no ",
    "fractional or decimal parts")

  ## Pass an 'R' with a fractional portion
  res <- try(bootES(data.frame(scores=1), R=0.5), silent=TRUE)
  expect_true(grepl(r_check_message, res))

  ## Pass an 'R' with string
  res <- try(bootES(data.frame(scores=1), R="5"), silent=TRUE)
  expect_true(grepl(r_check_message, res))

  ## Pass an 'R' with a negative value
  res <- try(bootES(data.frame(scores=1), R=-1), silent=TRUE)
  expect_true(grepl("R must be of length 1 and >= 0", res))
})

test_that("bootES handles invalid inputs", {
  gender = get_gender()

  g1 = c(11, 12, 13, 14, 15)
  g2 = c(26, 27, 28, 29)
  g3 = c(17, 18, 19, 20, 21, 22, 23)
  
  grpLabels = rep(c("A", "B", "C"), times=c(length(g1), length(g2), length(g3)))
  threeGps  = data.frame(grpLabels, scores=c(g1, g2, g3))
  threeGpsVec = c(g1, g2, g3)
  lambdas   = c(A=-1, B=2, C=-1)
    
  ## Pass a non-data.frame object as 'data'
  res = try(bootES("foo"), silent=TRUE)
  expect_true(grepl("'data' must be a data.frame or numeric vector.", res))

  ## Pass an invalid 'group' to 'contrast'
  res = try(bootES(gender, data.col="Meas3",
    group.col="Condition", contrast = c(Fake = -50, C = 50),
    scale.weights=TRUE), silent=TRUE)
  expect_true(grepl("'Fake' is/are not valid groups.", res[1]))
  
  ## Pass a data.frame to 'data' with no records
  expect_error(bootES(data.frame()))
  
  ## Pass an 'R' of length greater than 1
  expect_error(bootES(data.frame(scores=1), R=c(2, 1)))

  ## Use a 'data.col' not in 'data'
  expect_error(bootES(threeGps, R=250, data.col="foo"))

  ## Use a 'glass.control' value that is not a valid group
  res <- try(bootES(data.frame(scores=1), glass.control="foo"), silent=TRUE)
  expect_true(grepl("'glass.control' is not", res))
  
  ## Assert that using a block.col and glass.control together generates
  ## an error
  res <- try(bootES(gender, 
                    data.col="Meas1",
                    block.col="Gender", group.col="Condition",
                    glass.control="female",
                    contrast=c(A=1, B=-0.5, C=-0.5)), 
             silent=TRUE)
  expect_true(grepl("Cannot use 'block.col'", res))
  
  ## Use 'glass.control' and 'block.col' together w/out a data column
  ## and w/out contrasts
  res <- try(bootES(gender, block.col="Gender", grp.col="Condition"), 
             silent=TRUE)
  
  ## Use 'glass.control' and 'block.col' together w/out a data column
  res <- try(bootES(gender, block.col="Gender", group.col="Condition",
                    contrast=c(A=1, B=-0.5, C=-0.5)), 
             silent=TRUE)
  
  ## Assert that user cannot pass '...' arguments that are not valid 'boot'
  ## arguments.
  res <- try(bootES(gender, R=10, data.col="Meas1", block.col="Condition", 
                    slop.levels=letters[1:3]), 
             silent=TRUE)
  expect_true(grepl('invalid argument.*slop\\.levels', res[1], ignore.case=TRUE))

})

test_that("univariate statistics produce known results with bootES", {

  gender = get_gender()
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
  mean.res = bootES(threeGps, R=250, data.col="scores",
    effect.type="unstandardized")
  mean.res.vec = bootES(threeGpsVec, effect.type="unstandardized")
  expect_equal(truth, mean.res$t0)
  expect_equal(truth, mean.res.vec$t0)

  ## Test: 'rMeanBoot' through 'bootES'
  set.seed(1)
  truth     = bootES:::rMean(threeGps$scores)
  rMean.res = bootES(threeGps, R=300, data.col="scores", effect.type="r")
  expect_equal(truth, rMean.res$t0)

  ## Test: 'dMeanBoot' through 'bootES'
  set.seed(1)
  truth     = bootES:::dMean(threeGps$scores)
  rMean.res = bootES(threeGps, R=300, data.col="scores",
    effect.type="cohens.d")
  expect_equal(truth, rMean.res$t0)

  ## Test: 'dMeanBoot' and Cohen's Sigma d through 'bootES'
  set.seed(1)
  truth     = bootES:::dSigmaMeanBoot(threeGps$scores,
    1:length(threeGps$scores))
  rMean.res = bootES(threeGps, R=300, data.col="scores",
    effect.type="cohens.d.sigma")
  expect_equal(truth, rMean.res$t0)
  
  ## Test: 'hMeanBoot' and Hedge's g through 'bootES'
  set.seed(1)
  truth     = bootES:::hMean(threeGps$scores)
  rMean.res = bootES(threeGps, R=300, data.col="scores",
    effect.type="hedges.g")
  expect_equal(truth, rMean.res$t0)
  
  ## Test: 'akpRobustD' through 'bootES'
  dat = read.csv(system.file("robust_d_test.csv", package="bootES"))
  truth     = 0.190912
  set.seed(1)
  akp.res = bootES(dat, R=250, data.col="diff", effect.type="akp.robust.d")
  expect_equal(truth, akp.res$t0, tolerance=1e-4)
  
  set.seed(1)
  akp.res.1 = bootES(dat[["diff"]], R=250, effect.type="akp.robust.d")
  expect_equal(truth, akp.res.1$t0, tolerance=1e-4)
})

test_that("multivariate statistics produce known results with bootES", {

  gender = get_gender()
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  twoGpsA   = data.frame(x=c(g1, g2), team=grpLabels)
  twoGpsErr = data.frame(x=c(g1, g2), team=rep("A", length(c(g1, g2))))
  lambdas   = c(A=1, B=-1)
  
  ## Integration test of stat='contrast' and effect.type='unstandardized'
  set.seed(1)
  truth     = mean(g1) - mean(g2)
  unstdDiff.res = bootES(twoGpsA, R=250, data.col="x", group.col="team",
    effect.type="unstandardized", contrast=lambdas)
  expect_equal(truth, unstdDiff.res$t0)
   
  ## Integration test of stat='contrast' and effect.type='unstandardized' where
  ## there is only one group. This should cause an error.
  unstdDiff.err = try(bootES(twoGpsErr, R=250, data.col="x", group.col="team",
    effect.type="unstandardized"), silent=TRUE)    

  ## Integration test of stat='contrast' and effect.type='cohens.d.sigma'
  set.seed(1)
  test = bootES(gender, R=250, data.col="Meas1", group.col="Gender",
    effect.type="cohens.d.sigma", contrast=c(female=1, male=-1))
  expect_equal(-0.50104, test$t0, tolerance=1e-2)

  set.seed(1)
  test = bootES(gender, R=250, data.col="Meas1", group.col="Gender",
    effect.type="cohens.d.sigma", contrast=c(female=-1, male=+1))
  expect_equal(0.50104, test$t0, tolerance=1e-2)

  ## Integration test of stat='contrast' and effect.type='cohens.d.sigma'
  ## w/ glass control
  set.seed(1)
  test = bootES(gender, R=250, data.col="Meas1", group.col="Gender",
    effect.type="cohens.d.sigma", contrast=c('female', 'male'),
    glass.control='female')
  expect_equal(0.588, test$t0, tolerance=1e-2)

  ## Integration test of stat='contrast' and effect.type='cohens.d'
  ## w/ glass control
  set.seed(1)
  test = bootES(gender, R=250, data.col="Meas1", group.col="Gender",
    effect.type="cohens.d", contrast=c('female', 'male'),
    glass.control='female')
  expect_equal(0.557, test$t0, tolerance=1e-2)
  
  ## Integration test of stat='cor'
  set.seed(1)
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  twoGps    = data.frame(g1=g1)
  twoGps$g2 = rep(g2, length.out=nrow(twoGps))
  truth     = with(twoGps, cor(g1, g2))
  cor.res   = suppressWarnings(bootES(twoGps, R=10, effect.type="r"))
  expect_equal(truth, cor.res$t0)
})

test_that("hedges.g produces known results ", {
  gender = get_gender()

  set.seed(1)
  truth = 0.5096
  test  = bootES(gender, R=250, data.col = "Meas1", group.col = "Gender", 
                 contrast = c("female","male"), 
                 effect.type="hedges.g", glass.control="female")
  expect_equal(truth, test$t0, tolerance=1e-3)
})

test_that("'contrast' functionality produces known results", {
  gender = get_gender()

  ## Assert: Calculated value matches known value for an unstandardized
  ## contrast
  set.seed(1)
  truth = -522.43
  test  = bootES(gender, R=250, data.col="Meas3", group.col="Condition",
    contrast = c(A = -40, B = -10, C = 50), scale.weights=FALSE)
  expect_equal(truth, test$t0, tolerance=1e-2)

  ## Assert: Calculated value matches known value for an unstandardized
  ## contrast with weights scaled
  set.seed(1)
  truth.contrast.scaled = -10.4486
  test  = bootES(gender, R=250, data.col="Meas3", group.col="Condition",
    contrast = c(A = -40, B = -10, C = 50), scale.weights=TRUE)
  expect_equal(truth.contrast.scaled, test$t0, tolerance=1e-4)

  ## Assert: Calculated value matches known value for an unstandardized
  ## contrast with weights scaled and a group left out
  set.seed(1)
  truth.contrast.omit = -3.0535
  test  = bootES(gender, R=250, data.col="Meas3", group.col="Condition",
    contrast = c(A = -1, C = 1))  
  expect_equal(truth.contrast.omit, test$t0, tolerance=1e-4)

  ## Assert: Default weights of -1 and 1 are used when not passed in
  test.dflt  = bootES(gender, R=250, data.col="Meas3", group.col="Condition",
    contrast = c('A', 'C'))
  expect_equal(truth.contrast.omit, test.dflt$t0, tolerance=1e-4)

  ## Assert: Scales user-specified weights
  set.seed(1)
  truth = -0.45488
  test = bootES(gender, R=250, data.col = "Meas1", group.col = "Gender",
    contrast=c(female = 3, male = -3), effect.type = "hedges.g")
  expect_equal(truth, test$t0, tolerance=1e-3)
  
  gender$GenderByCond = paste(gender$Gender, gender$Condition, sep = "-")
  set.seed(1)
  truth = 46.71499
  test <- bootES(gender, data.col="Meas1", group.col="GenderByCond", 
                 contrast = c("female-A" = -40, "male-A" = -40, 
                              "female-B" = -10, "male-B" = -10, 
                              "female-C" = 50, "male-C" = 50))
  expect_equal(truth, test$t0, tolerance=1e-3)

  ## Assert: Test the blocking column w/ contrasts that must be scaled
  set.seed(1)
  truth = 46.71499
  test = bootES(gender, R=250, data.col="Meas1",
    block.col="GenderByCond", group.col="Condition",
    contrast=c(A=-40, B=-10, C=50))
  expect_equal(truth, test$t0, tolerance=1e-3)
  
  ## Assert: Test the blocking column w/ contrasts that don't need to
  ## be scaled
  set.seed(1)
  truth = 36.783
  test <- bootES(gender, R=250, data.col="Meas1", group.col="Gender",
                 contrast=c("female"=-1, "male"=1), block.col="Condition")
  expect_equal(truth, test$t0, tolerance=1e-3)
  cond_means <- with(gender, tapply(Meas1, GenderByCond, mean))
  
  ## Assert: Test that means are unweighted
  set.seed(1005)
  # means <- with(gender, tapply(Meas1, Condition, mean))
  # truth <- with(gender, mean(tapply(Meas1, Condition, mean)))
  truth <- 266.944 # unweighted mean 
  test <- bootES(gender, R=250, data.col="Meas1", block.col="Condition")
  expect_equal(truth, test$t0, tolerance=1e-3)
#   wt_mean <- local({
#     means <- with(gender, tapply(Meas1, GenderByCond, mean))
#     counts <- with(gender, tapply(Meas1, GenderByCond, length))
#     sum(means * counts / nrow(gender))
#   })

})

test_that("apk.robust.d produces known results", {
  gender = get_gender()
  truth = 0.5487
  test.1 = bootES(gender, R=250, data.col="Meas1", group.col="Gender",
                  contrast=c(male=1, female=-1), 
                  effect.type="akp.robust.d")
  expect_equal(test.1$t0, truth, tolerance=1e-3)
})

test_that("'cor.dirr' functionality produces known results", {
  ## Unit test of stat='cor.diff'
  ## Note that R=300 set so runs quickly but does not generate warnings
  ## when calculating CIs
  gender = get_gender()
  set.seed(1)

  iris = get(data("iris"))
  cols = c("Species", "Sepal.Length", "Petal.Length")
  iris_sv = iris[iris$Species %in% c("setosa", "versicolor"), cols]
  iris_ls = split(iris_sv, iris_sv$Species)
  setosa = iris_ls[["setosa"]]
  versicolor = iris_ls[["versicolor"]]
  truth = (cor(setosa$Sepal.Length, setosa$Petal.Length) -
    cor(versicolor$Sepal.Length, versicolor$Petal.Length))

  cor.diff.res = bootES(iris_sv, R=300, group.col="Species", effect.type="r")
  expect_equal(truth, cor.diff.res$t0)
})

test_that("slope functionality produces known results", {
  ## Regression test for when effect.type='slope'
  gender = get_gender()
  set.seed(1)
  truth <- -0.1244
  test  <- bootES(gender, R=200, data.col="Meas3", group.col="Condition",
                  slope.levels=c(A=30, B=60, C=120))
  expect_equal(truth, test$t0, tolerance=1e-2)
  
  test  <- bootES(gender, R=200, data.col="Meas3", slope.levels="Dosage")
  expect_equal(truth, test$t0, tolerance=1e-2)
                  
})

test_that("'effect.type' functionality produces known results", {
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
  mean.res = bootES(threeGps, R=250, data.col="scores",
    effect.type="unstandardized")
  mean.res.vec = bootES(threeGpsVec, effect.type="unstandardized")
  expect_equal(truth, mean.res$t0)
  expect_equal(truth, mean.res.vec$t0)

  ci.types = eval(formals(bootES)$ci.type)
  ci.types = ci.types[!ci.types %in% "stud"]
  for (ci.type in ci.types) {
    set.seed(1)
    if (ci.type == "stud") {
      . <- bootES(threeGps, R=250, data.col="scores",
                  effect.type="unstandardized", ci.type=ci.type,
                  var.t0=1, var.t=1)
    } else {
      . <- bootES(threeGps, R=250, data.col="scores",
                  effect.type="unstandardized", ci.type=ci.type)
    }
  }
})

test_that("blocking functionality produces known results", {
  ## testgroup: Calls calcUnstandardizedMean through the bootES interface
  
  ## Assert that blocking and grouping work exactly the same when the 
  ## contrasts are specified at the block or group level

  gender = get_gender()

  set.seed(1)
  test.1a = bootES(gender, R=999, data.col="Meas1", group.col="Gender", 
                   contrast=c(female=-1, male=1), block.col="GenderByCond")
  
  set.seed(1)
  test.1b = bootES(gender, R=999, 
                   data.col = "Meas1", group.col = "GenderByCond", 
                   contrast = c("female-A"=-1, "female-B"=-1, "female-C"=-1, 
                                "male-A"=1, "male-B"=1, "male-C"=1))
  
  expect_equal(summary(test.1a), summary(test.1b))
  
  ## testgroup: Calls calcCohensD through the bootES interface
  set.seed(1)
  test.2a = bootES(gender, R=999, data.col="Meas1", group.col="Gender", 
                   contrast=c(female=-1, male=1), block.col="GenderByCond",
                   effect.type="cohens.d")
  
  set.seed(1)
  test.2b = bootES(gender, R=999, effect.type="cohens.d",
                   data.col = "Meas1", group.col = "GenderByCond", 
                   contrast = c("female-A"=-1, "female-B"=-1, "female-C"=-1, 
                                "male-A"=1, "male-B"=1, "male-C"=1))
  expect_equal(summary(test.1a), summary(test.1b))
  
  ## test: Assert that bootES automatically crosses them
  set.seed(1)
  test.3 = bootES(gender, R=999, data.col="Meas1", group.col="Gender", 
                  contrast=c(female=-1, male=1), block.col="Condition")
  expect_equal(summary(test.3), summary(test.1a))
  
})
