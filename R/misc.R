## Daniel Gerlanc and Kris Kirby (2010)
## Miscellaneous helper functions

calcSlopeLambdas <- function(vals, grps, contrasts) {
  ## Calculates the weights to use for a slope calculation
  ##
  ## Args:
  ##   vals: a numeric vector
  ##   grps: a numeric, character, or factor vector of groups
  ##   contrasts: a named, numeric vector mapping 'grps' to predictor values
  ##              otherwise NULL
  ##
  ## Returns:
  ##  a vector of weights
  
  avgs = c(tapply(vals, grps, mean))
  preds <- if (is.numeric(contrasts))
    contrasts[names(avgs)] else as.numeric(names(avgs))
  ss    = sum((preds - avgs)^2)
  lmbds = (preds - avgs) / ss
}

scaleLambdasBySide <- function(lambdas) {
  ## Scales negative and positive lambdas to sum to +/- 1, respectively
  ##
  ## Args:
  ##   lambdas: A numeric vector of lambdas
  ##
  ## Returns:
  ##  a scaled numeric vector of lambdas

  signs = sign(lambdas)
  sum.by.side = c(tapply(lambdas, signs, sum))
  neg.idx = signs == -1 
  pos.idx = signs == +1
  
  lambdas[neg.idx]  = lambdas[neg.idx] / abs(sum.by.side["-1"])
  lambdas[pos.idx]  = lambdas[pos.idx] / sum.by.side["1"]
  return(lambdas)
}

pooledSD <- function(vals, grps, pop.sd=FALSE) {
  ## Calculate the pooled standard deviation of a vector.
  ## Args:
  ##  @param x a numeric vector
  ##  @param grps a numeric, character, or factor vector specifying groups
  ##  @param pop.sd calculate the biased population standard deviation instead
  ##         of the unbiased sample standard deviation
  ##
  ##  @return
  ##    An unbiased point estimate of the pooled standard deviation

  grps.ls = split(vals, grps, drop=TRUE)
  ss      = sum(sapply(grps.ls, function(x) sum((x - mean(x))^2)))
  denom   = if (pop.sd) length(vals) else length(vals) - length(grps.ls)
  res     = sqrt(ss / denom)
  return(res)
}
