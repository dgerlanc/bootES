## Daniel Gerlanc and Kris Kirby (2010)
## Input and regression tests for the scaleLambdasBySide function

test.scaleLambdasBySide <- function() {
  ## Test the functioning of the scaleLambdasBySide function

  ## Simplest Case
  x1   = c(-1, 1)
  res1 = bootES:::scaleLambdasBySide(x1)
  checkEquals(x1, res1)

  x2   = c(-0.33, 0.33)
  res2 = bootES:::scaleLambdasBySide(x2)
  checkEquals(res2, c(-1, 1))

  x3   = c(-0.25, -0.5, 0.5)
  res3 = bootES:::scaleLambdasBySide(x3)
  truth3 = c(-1/3, -2/3, 1)
  checkEquals(res3, truth3, tol=1e-8)
  
  x4     = c(-5, 0, 6)
  res4   = bootES:::scaleLambdasBySide(x4)
  truth4 = c(-1, 0, 1)
  checkEquals(res4, truth4)

  x5     = c(0, 0)
  res5   = bootES:::scaleLambdasBySide(x5)
  checkEquals(res5, x5)
  
}
