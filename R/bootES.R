## Daniel Gerlanc and Kris Kirby (2010)
## High-level function for bootstrap analyses using the 'boot' package

bootES <- function(dat, R=1000, data.col, grp.col,
                   stat=c("mean", "contrast", "cor", "cor.diff", "slope"),
                   effect.type=c("unstandardized", "cohens.d", "hedges.g",
                     "cohens.d.sigma", "r"),
                   contrasts=NULL,                   
                   glass.control=NULL,
                   scale.weights=FALSE,
                   ci.type=c("bca", "norm", "basic", "stud", "perc", "all", "none"),
                   ci.conf=0.95,
                   verbose=0) {
  
  ## Performs different variants of bootstrap analyses for calculating effect sizes.
  ##
  ## Args:
  ##   dat           : a vector or data frame containing the one or more columns of values
  ##                   (required), group labels (optional), and contrasts (optional)
  ##                   for each sample
  ##   R             : the number of bootstrap 'repetitions' to perform
  ##   data.col      : The column in 'dat' containing the sample values
  ##   grp.col       : The column in 'dat' containing the grouping info
  ##   stat          : The statistic to calculate
  ##   effect.type   : The type of standardization to perform when calculating 'stat'
  ##   contrasts     : A named vector specifying the lambdas for different groups in 'dat'.
  ##                   The default applies to the simplest, 2-group unweighted case
  ##   glass.control : The group for which the standard deviation should be used, eg.
  ##     "glass.control='A'"
  ##   scale.weights : TRUE/FALSE, scale the lambdas to [-1, 1]
  ##   ci.type       : The type of confidence interval to generate (see 'boot.ci')
  ##   ci.conf       : The confidence level of the interval
  ##   verbose       : Higher levels generate more output
  ##
  ## Returns:
  ##   An object of class 'bootES' and 'boot'
  ##
  ## Details: If 'R' is not a whole number, it will be round down to the nearest whole number.
  ##   * stat=cor: 'dat' must be a two-column data frame, where each of the columns is numeric
  ##

  ## Error handling
  stat        = match.arg(stat)
  effect.type = match.arg(effect.type)
  ci.type     = match.arg(ci.type)
  
  ## Checks on 'dat'.
  if (!(is.data.frame(dat) || is.numeric(dat)))
    stop("'dat' must be a data.frame or numeric vector.")

  if (is.numeric(dat)) {
      dat      = data.frame(scores=dat, row.names=NULL)
      data.col = "scores"
  }
  
  if (!nrow(dat) > 0)
    stop("'dat' contains no records!")
  
  ## Checks on 'R'.
  R = as.integer(R)
  r.is.valid = (length(R) == 1) && is.numeric(R) && R > 0
  if (!r.is.valid)
    stop("R must be an integer of length 1 and greater than 0")
  
  ## Check and extract 'data.col'.
  if (!missing(data.col)) {
    if (!is.character(data.col))
      stop("'data.col' must be a character vector.")
  
    if (!data.col %in% colnames(dat))
      stop("'data.col' missing from 'dat'")

    vals = dat[[data.col]]
  }
  
  ## Check and extract 'grp.col'.
  grps = NULL
  if (!missing(grp.col)) {
    if (!is.character(grp.col)) 
      stop("'grp.col' must be a character vector.")
    
    if (!grp.col %in% colnames(dat))
      stop("'grp.col' missing from 'dat'")
    
    grps = as.factor(dat[[grp.col]])
  }  

  ## Checks on scale.weights
  if (!is.logical(scale.weights) || length(scale.weights) != 1)
    stop("'scale.weights' must be a logical vector of length 1.")
  
  ## Check and extract contrasts.
  lmbds = NULL
  if (!missing(contrasts)) {
    if (is.null(names(contrasts)))
      stop("'contrasts' must be a named vector")
    
    if (missing(grp.col))
      stop("Must specify a 'grp.col' when providing a 'contrasts' argument.")

    ## Scale contrasts if specified and not using the slope function
    if (scale.weights && abs(sum(contrasts)) > 1e-4 && stat != "slope")
      lmbds = scaleLambdasBySide(lmbds)    
    
    lmbds = contrasts[match(grps, names(contrasts))]
    
    ## Check that there are no NA lambdas.
    na.lmbds = is.na(lmbds)
    no.lmbds = sort(unique(grps[na.lmbds]))
    if (any(na.lmbds))
      stop(paste(sQuote(no.lmbds), collapse=", "), " have missing lambdas.")
  }
  
  ## Error handling for the stat='cor'
  if (stat == "cor") {
    is.valid = length(dat) == 2 && all(sapply(dat, is.numeric))
    if (!is.valid)
      stop("'dat' must be a data frame with two numeric columns.")
  }

  ## Error handling for the stat='cor.diff'
  if (stat == "cor.diff") {
    is.valid = is.numeric(dat[[1]]) && is.numeric(dat[[2]])
    if (!is.valid)
      stop("'The first two columns of 'dat' must be numeric.'")

    if (missing(grp.col))
      stop("You must specify a grouping column.")
    
    if (length(unique(grps)) != 2)
      stop("There must be only 2 distinct groups!")
  }
  
  ## Error handling for slope argument
  if (stat == "slope") {
    if (is.null(vals))
      stop("Invalid data column.")
      
    if (is.null(grps))
      stop("Invalid grouping column.")

    if (!is.numeric(grps)) {
      err.msg <- paste("If data in group column is not numeric, a numeric",
                       "'contrasts' argument with names corresponding to the",
                       "values in the 'grp.col' must be provided.", collapse=" ")
      if (!is.numeric(contrasts)) 
        stop(err.msg)

      if (!all(names(contrasts) %in% grps))
        stop(err.msg)
    } else {
      contrasts = NULL
    }
    
    if (length(unique(grps)) < 2)
      stop("There must be at least 2 distinct groups!")
  }
  
  ## Error handling for 'glass.control'
  if (!missing(glass.control)) {
    if (!is.character(glass.control) || length(glass.control) != 1)
      stop("glass.control must be a character vector of length 1.")

    if (!glass.control %in% colnames(dat))
      stop("'glass.control' missing from 'dat'")
  }  

  ## Error handling for 'verbose'
  verbose = as.integer(verbose)
  if (!verbose %in% 0:2)
    stop("'verbose' must be 0, 1, or 2.")
  
  ## Simplest Case: No groups, so we can calculate all of the stats for a single group
  n.grps = length(unique(grps))
  single.group = (is.null(grps) || n.grps == 1) && stat != "cor"
  if (single.group) {
    boot.fun <- switch(effect.type,
                       unstandardized = meanBoot,
                       r              = rMeanBoot,
                       hedges.g       = dMeanBoot,
                       cohens.d       = dMeanBoot,
                       cohens.d.sigma = dSigmaMeanBoot)
    res = boot(vals, statistic=boot.fun, R=R) # TODO: Parse, then evaluate this call
  } else { # Two or more groups
    if (effect.type %in% c("cohens.d", "hedges.g", "cohens.d.sigma")) {
      res = boot(vals, calcCohensD, R, stype="f", strata=grps, grps=grps,
        lambdas=lmbds, cohens.d.sigma=(effect.type == "cohens.d.sigma"),
        glass.control=glass.control)
    } else if (stat == "cor") {
      res = boot(dat, statistic=corBoot, R=R) 
    } else if (stat == "cor.diff") {
      res = boot(dat, calcBootCorDiff, R=R, stype="f", strata=grps, grps=grps)
    } else if (stat == "slope") {      
      lmbds = calcSlopeLambdas(vals, grps, contrasts)
      res   = boot(vals, calcUnstandardizedMean, R=R, stype="f", strata=grps,
        grps=grps, lambdas=lmbds)
    } else if (stat == "mean" && effect.type == "unstandardized") {
      res = boot(vals, meanUnweightedBoot, R, stype="f", strata=grps, grps=grps)
    } else if (stat == "contrast" && effect.type == "unstandardized") {
      res = boot(vals, meanDiffBoot, R, stype="f", strata=grps, grps=grps)
    } else if (stat == "contrast" || effect.type == "r") {
      res = boot(vals, calcPearsonsR, R, stype="f", strata=grps, grps=grps,
        lambdas=lmbds)
    } else {
      stop("This combination of 'stat' and 'effect.type' not yet implemented!")
    }
  }

  res[["ci.type"]] = ci.type
  res[["ci.conf"]] = ci.conf
  res[["verbose"]] = verbose
  class(res) = c("bootES", "boot")
  
  if (verbose > 0L)
    printTerse(res)

  if (verbose > 1L)
    plot(res)
  
  return(res)
}

print.bootES <- function(x, ...) {
    ## Prints the bootstrap summary statistics for a bootES object, followed by
    ## the confidence interval if the user required its calculation.
    class(x) <- "boot"            
    boot:::print.boot(x)
    cat("\n")
    if (x[["ci.type"]] != "none")
      print(boot.ci(x, conf=x[["ci.conf"]], type=x[["ci.type"]]))
}

printTerse <- function(x) {
  ## Print a terse representation of the bootES results describing the
  ## confidence level, the type of confidence interval, the statistic, and the
  ## calculated CI bounds
  stopifnot(inherits(x, "bootES"))
  
  ## Extract the confidence interval
  ci.type = x[["ci.type"]]
  stat = x[["t0"]]
  ci = boot.ci(x, conf=x[["ci.conf"]], type=ci.type)
  ci = ci[[ci.type, exact=FALSE]]

  bounds = switch(ci.type,
    bca   = ci[1, 4:5, drop=TRUE],
    norm  = ci[1, 2:3, drop=TRUE],
    all   = ci[["bca"]][1, 4:5, drop=TRUE],
    ci)

  nms = c("Stat", "CI (Low)", "CI (High)")
  res = matrix(c(stat, bounds), nrow=1, dimnames=list(NULL, nms))

  cat(sprintf("%.2f%% %s Confidence Interval\n", 100 * x[["ci.conf"]], ci.type))
  print(res)
  cat("\n")
}

determineStat <- function(dat, grps=NULL, effect.type=NULL, contrasts=NULL) {
  ## Based on the arguments passed to bootES, determine the type of statistic to
  ## calculate
  res = ''
  if (is.null(grps)) {
    if (ncol(dat) == 1) 
      res = 'mean'
    else
      res <- if (identical(effect.type, 'slope')) 'slope' else 'cor'
  } else {
    if (!is.null(contrasts)) {
      res = 'contrast'
    } else {
      if (length(unique(grps)) == 2)
        res = 'cor.diff'
      else
        stop('Within-subject difference between correlations not implemented.')
             
    }
  }
  return(res)
}
