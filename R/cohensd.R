## Daniel Gerlanc and Kris Kirby (2010)
## Functions for calculating Cohen's D

calcCohensD <- function(vals, freq, grps,
                        cohens.d.sigma=FALSE,
                        glass.control="") {
  ## Compute Cohen's d/Hedge's g or optionally Cohen's Sigma D for the
  ## mean of one or more groups
  ## 
  ## Args:
  ##   @param vals a numeric vector of values
  ##   @param freq a frequency vector of length equal to the number of records
  ##     in 'dat' indicating how many times an observation should be drawn
  ##     from each sample.
  ##   @param grps a grouping vector of the same length as 'vals'
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

    if (identical(glass.control, nm))
      glass.sd = sqrt(ss[nm] / (length(sample.vals) - 1))
  }

  ## Calculate the standard deviation.
  if (!is.null(glass.sd)) {
    sd.hat = glass.sd
  } else {
    num.grps     = length(grp.idx)
    sd.hat.denom = if (cohens.d.sigma) length(vals) else length(vals) - num.grps
    sd.hat       = sqrt(sum(ss) / sd.hat.denom)
  }
  
  res = (means[1] - means[2]) / sd.hat
  res = as.vector(res) # remove names from 'means'
  return(res)
}
