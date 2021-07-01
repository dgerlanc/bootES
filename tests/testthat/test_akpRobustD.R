## Daniel Gerlanc and Kris Kirby (2010-2015)
## Input and regression tests for the akpRobustD function

test_that("akpRobustD produces a known result", {
  dat = read.csv(system.file("robust_d_test.csv", package="bootES"))
  truth = 0.190912
  res = bootES:::akpRobustD(dat[["diff"]], seq_len(nrow(dat)))
  expect_equal(truth, res, tolerance=1e-4)
})
