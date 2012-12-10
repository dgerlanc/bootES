## Daniel Gerlanc and Kris Kirby (2010-2012)
## Functions for calculating Cohen's D

calcCohensD <- function(vals, freq, grps,
                        contrast,
                        hedges.g=FALSE,
                        cohens.d.sigma=FALSE,
                        glass.control="") {
  ## Compute Cohen's d/Hedge's g or optionally Cohen's Sigma D for the
  ## mean of one or more groups
  ## 
  ## Args:
  ##   @param vals a numeric vector of values
  ##   @param freq a frequency vector of length equal to the number of records
  ##     in 'data' indicating how many times an observation should be drawn
  ##     from each sample.
  ##   @param grps a grouping vector of the same length as 'vals'
  ##   @param contrast a named, numeric vector of contrast weights.
  ##     The names must match the values in 'grps'
  ##   @param hedges.g Whether to make the hedges.g adjustment to cohens.d
  ##   @param cohens.d.sigma Whether to use the population standard
  ##     deviation instead of the sample standard deviation
  ##   @param glass.control a character vector of length 1 specifying the
  ##     optional group to use to calculate standard deviation
  ##
  ## Details:
  ##   This function is meant to be passed as the 'statistic' argument to the
  ##   'boot' function. 'freq' should be a frequency vector of the type returned
  ##   by 'boot'.

  ## Get the integer indices of the different groups.
  grp.idx = split(seq_along(vals), grps, drop=TRUE)   
  grp.nms = names(grp.idx)
  
  ## Calculate the means and sums-of-squares for the bootstrap samples from each
  ## group.
  means = ss = numeric()
  glass.sd = NULL
  for (nm in grp.nms) {
    this.grp.idx = grp.idx[[nm]]
    sample.vals  = rep(vals[this.grp.idx], times=freq[this.grp.idx])
    means[nm]    = mean(sample.vals)
    ss[nm]       = sum((sample.vals - means[nm])^2)

    if (identical(glass.control, nm)) {
      if (cohens.d.sigma) {
        denom = length(sample.vals)
      } else {
        denom = length(sample.vals) - 1
      }
      glass.sd = sqrt(ss[nm] / denom)
    }
  }

  ## Calculate the standard deviation.
  num.grps     = length(grp.idx)
  df = if (cohens.d.sigma) length(vals) else length(vals) - num.grps
  if (!is.null(glass.sd)) {
    sd.hat = glass.sd
  } else {    
    sd.hat = sqrt(sum(ss) / df)
  }

  res = sum(contrast[grp.nms] * means) / sd.hat
  res = as.vector(res) # remove names from 'means'
  if (hedges.g) {
    hadj = gamma(df/2)/( sqrt(df/2)*gamma((df-1)/2) )
    res  = hadj * res
  }
  return(res)
}
