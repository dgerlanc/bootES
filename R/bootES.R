## Daniel Gerlanc and Kris Kirby (2010)
## High-level function for bootstrap analyses using the 'boot' package

bootES <- function(dat, R=1000, data.col=NULL, grp.col=NULL,
                   effect.type=c("unstandardized", "cohens.d", "hedges.g",
                     "cohens.d.sigma", "r", "slope"),
                   contrasts=NULL,                   
                   glass.control=NULL,
                   scale.weights=TRUE,
                   ci.type=c("bca", "norm", "basic", "stud", "perc",
                     "all", "none"),
                   ci.conf=0.95,
                   verbose=0) {
  
  ## Performs different variants of bootstrap analyses for calculating
  ## effect sizes.
  ##
  ## Args:
  ##   dat           : a vector or data frame containing the one or more
  ##                   columns of values (required), group labels (optional),
  ##                   and contrasts (optional) for each sample
  ##   R             : the number of bootstrap 'repetitions' to perform
  ##   data.col      : The column in 'dat' containing the sample values
  ##   grp.col       : The column in 'dat' containing the grouping info
  ##   effect.type   : The type of standardization to perform
  ##   contrasts     : A named vector specifying the lambdas for different
  ##                   groups in 'dat'
  ##   glass.control : The group for which the standard deviation should
  ##                   be used, eg. "glass.control='A'"
  ##   scale.weights : TRUE/FALSE, scale the lambdas to [-1, 1]
  ##   ci.type       : The type of confidence interval to generate
  ##                   (see 'boot.ci')
  ##   ci.conf       : The confidence level of the interval
  ##   verbose       : Higher levels generate more output
  ##
  ## Returns:
  ##   An object of class 'bootES' and 'boot'
  ##
  ## Details: If 'R' is not a whole number, it will be round down to the nearest
  ## whole number.  * stat=cor: 'dat' must be a two-column data frame, where
  ## each of the columns is numeric
  ##
  ## TODO: Parse, then evaluate the 'boot' call

  ## Error handling
  effect.type = match.arg(effect.type)
  ci.type     = match.arg(ci.type)
  
  ## Checks on 'dat'.
  if (!(is.data.frame(dat) || is.numeric(dat)))
    stop("'dat' must be a data.frame or numeric vector.")

  ## If 'dat' has been passed as a numeric vector, save it as a data.frame
  ## with a data.col='scores'
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
  if (!is.null(data.col)) {
    if (!is.character(data.col))
      stop("'data.col' must be a character vector.")
  
    if (!data.col %in% colnames(dat))
      stop("'data.col' missing from 'dat'")

    vals = dat[[data.col]]
  }
  
  ## Check and extract 'grp.col'.
  grps = NULL
  if (!is.null(grp.col)) {
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
  lmbds.orig = contrasts
  if (!is.null(contrasts)) {
    if (is.null(names(contrasts)))
      stop("'contrasts' must be a named vector")
    
    if (is.null(grp.col))
      stop("Must specify a 'grp.col' when providing a 'contrasts' argument.")

    lmbds = contrasts
    
    ## Scale contrasts if specified and not using the slope function
    scale.lambdas = scale.weights && effect.type != "slope"
    if (scale.lambdas)
      lmbds = scaleLambdasBySide(lmbds)    
    
    ## Assert that there are no NA lambdas, then subset 'dat' to the groups for
    ## which contrasts were provided.
    lmbds.exist   = names(lmbds) %in% grps
    missing.lmbds = names(lmbds)[!lmbds.exist]
        
    if (length(missing.lmbds))
      stop(paste("'", missing.lmbds, "'", sep="", collapse=", "),
           " is/are not valid groups.")
    
    boot.groups = unique(names(lmbds))
    dat  = dat[dat[[grp.col]] %in% boot.groups, ]

    vals = dat[[data.col]]
    grps = as.factor(dat[[grp.col]])
  }

  ## Determine the 'stat' based on the passed arguments
  stat = determineStat(dat, data.col, grps, effect.type, lmbds)
  
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

    if (is.null(grp.col))
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
                       "values in the 'grp.col' must be provided.",
                       collapse=" ")
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
  if (!is.null(glass.control)) {
    if (!is.character(glass.control) || length(glass.control) != 1)
      stop("glass.control must be a character vector of length 1.")

    if (!glass.control %in% colnames(dat))
      stop("'glass.control' missing from 'dat'")
  }  

  ## Error handling for 'verbose'
  verbose = as.integer(verbose)
  if (!verbose %in% 0:2)
    stop("'verbose' must be 0, 1, or 2.")
  
  ## Simplest Case: No groups, so we can calculate all of the stats for a single
  ## group
  n.grps = length(unique(grps))
  single.group = (is.null(grps) || n.grps == 1) && stat != "cor"
  if (single.group) {
    boot.fun <- switch(effect.type,
                       unstandardized = meanBoot,
                       r              = rMeanBoot,
                       hedges.g       = dMeanBoot,
                       cohens.d       = dMeanBoot,
                       cohens.d.sigma = dSigmaMeanBoot)    
    res = boot(vals, statistic=boot.fun, R=R) 
  } else { # Two or more groups
    if (effect.type %in% c("cohens.d", "hedges.g", "cohens.d.sigma")) {
      res = boot(vals, calcCohensD, R, stype="f", strata=grps, grps=grps,
        cohens.d.sigma=(effect.type == "cohens.d.sigma"),
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
      res = boot(vals, calcUnstandardizedMean, R, stype="f", strata=grps,
        grps=grps, lambdas=lmbds)
    } else if (stat == "contrast" || effect.type == "r") {
      res = boot(vals, calcPearsonsR, R, stype="f", strata=grps, grps=grps,
        lambdas=lmbds)
    } else {
      fmt = paste("effect.type: %s for a 'dat' of class '%s' of length %d",
        "with data.col of class %s and grp.col of class %s not implemented.",
        collapse="")
      msg = sprintf(msg, effect.type, class(dat), length(dat), class(grps),
        class(vals))
      stop(msg)
    }
  }

  res[["ci.type"]] = ci.type
  res[["ci.conf"]] = ci.conf
  res[["verbose"]] = verbose
  res[["contrasts"]] = lmbds.orig
  res[["contrasts.scaled"]] = lmbds
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

  ## Print the scaled and unscaled contrasts.
  if (!is.null(x$contrasts)) {
    cat(sprintf("User-specified lambdas: (%s)\n",
                paste(x$contrasts, collapse=", ")))
  }
  if (!is.null(x$contrasts.scaled)) {
    cat(sprintf("Scaled lambdas: (%s)\n",
                paste(x$contrasts.scaled, collapse=", ")))
  }
  
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

determineStat <- function(dat,
                          data.col=NULL,
                          grps=NULL,
                          effect.type=NULL,
                          contrasts=NULL) {
  ## Based on the arguments passed to bootES, determine the type of statistic to
  ## calculate
  ##
  ## Args:
  ##  dat:
  ##  grps:
  ##  effect.type:
  ##  contrasts:
  ## 

  do.cor = (is.data.frame(dat) && ncol(dat) == 2 && is.numeric(dat[[1]]) && 
            is.numeric(dat[[2]]))
  res = ''
  if (is.null(grps)) {
    ## Single Group
    ## - slope -> slope
    ## - data.col == NULL -> cor
    ## - r -> r
    ## - otherwise -> mean
    if (identical(effect.type, 'slope')) {
      res = 'slope'
    } else if (is.null(data.col) && do.cor) {
      res = 'cor'
    } else if (identical(effect.type, 'r')) {
      res = 'r'
    } else {
      res = 'mean'
    }
  } else { 
    if (!is.null(contrasts)) {
      res = 'contrast'
    } else {
      if (is.numeric(dat)) {
        res = 'mean'
      } else {
        useCorDiff = is.data.frame(dat) && length(unique(grps)) == 2 &&
        ncol(dat) > 2
        if (useCorDiff)
          res = 'cor.diff'
        else
          stop('Within-subject difference between correlations not implemented.')
      }
    }
  }
  return(res)
}
