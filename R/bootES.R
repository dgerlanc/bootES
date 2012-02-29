## Daniel Gerlanc and Kris Kirby (2010)
## High-level function for bootstrap analyses using the 'boot' package

bootES <- function(dat, R=1000, data.col=NULL, group.col=NULL,
                   effect.type=c("unstandardized", "cohens.d", "hedges.g",
                     "cohens.d.sigma", "r", "slope"),
                   contrasts=NULL,
                   slope.levels=NULL,
                   glass.control=NULL,
                   scale.weights=TRUE,
                   ci.type=c("bca", "norm", "basic", "perc", "stud", "none"),
                   ci.conf=0.95,
                   plot=FALSE,
                   ...) {

  
  ## Performs different variants of bootstrap analyses for calculating
  ## effect sizes.
  ##
  ## Args:
  ##   dat           : a vector or data frame containing the one or more
  ##                   columns of values (required), group labels (optional),
  ##                   and contrasts (optional) for each sample
  ##   R             : the number of bootstrap 'repetitions' to perform
  ##   data.col      : The column in 'dat' containing the sample values
  ##   group.col       : The column in 'dat' containing the grouping info
  ##   effect.type   : The type of standardization to perform
  ##   contrasts     : A named vector specifying the lambdas for different
  ##                   groups in 'dat'
  ##   slope.levels  : A named vector specifying the levels for different
  ##                   groups in 'dat'
  ##   glass.control : The group for which the standard deviation should
  ##                   be used, eg. "glass.control='A'"
  ##   scale.weights : TRUE/FALSE, scale the lambdas to [-1, 1]
  ##   ci.type       : The type of confidence interval to generate
  ##                   (see 'boot.ci')
  ##   ci.conf       : The confidence level of the interval
  ##   plot          : Generate the plot?
  ##   ...           : additional arguments passed to 'boot.ci'
  ##
  ## Returns:
  ##   An object of class 'bootES' and 'boot'
  ##
  ## Details: If 'R' is not a whole number, it will be round down to the nearest
  ## whole number.  * stat=cor: 'dat' must be a two-column data frame, where
  ## each of the columns is numeric

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
  
  ## Check and extract 'group.col'.
  grps = NULL
  if (!is.null(group.col)) {
    if (!is.character(group.col)) 
      stop("'group.col' must be a character vector.")
    
    if (!group.col %in% colnames(dat))
      stop("'group.col' missing from 'dat'")
    
    grps = as.factor(dat[[group.col]])
  }  

  ## Checks on scale.weights
  if (!is.logical(scale.weights) || length(scale.weights) != 1)
    stop("'scale.weights' must be a logical vector of length 1.")

  ## Process arguments for slope calculations
  lmbds = NULL
  lmbds.orig = contrasts
  if (effect.type == "slope") {
    if (is.null(vals))
      stop("Invalid 'data.col'.")
    
    if (is.null(grps))
      stop("Invalid 'group.col'")
    
    if (is.null(slope.levels))
      stop("Must specify 'slope.levels'.")

    if (!is.null(contrasts))
      stop("Cannot specify 'contrasts' and 'slope.levels'")
    
    invalid.levels <- !is.numeric(slope.levels) || is.null(names(slope.levels))
    if (invalid.levels)
      stop("'slope.levels' must be a named numeric vector.")
    
    lmbds <- calcSlopeLambdas(slope.levels)

    ## Assert that there are no NA slopes, then subset 'dat' to the groups for
    ## which slope.levels were provided.
    lmbds.exist   = names(lmbds) %in% grps
    missing.lmbds = names(lmbds)[!lmbds.exist]
        
    if (length(missing.lmbds))
      stop(paste("'", missing.lmbds, "'", sep="", collapse=", "),
           " is/are not valid groups.")
    
    boot.groups = unique(names(lmbds))
    dat  = dat[dat[[group.col]] %in% boot.groups, ]

    grps = as.factor(dat[[group.col]])
    vals = dat[[data.col]]
  }

  
  ## Check and extract contrasts.  
  if (!is.null(contrasts)) {            
    use.default.contrasts = is.character(contrasts) && length(contrasts) == 2
    if (use.default.contrasts)
      contrasts = structure(c(-1, 1), names=contrasts)

    if (is.null(names(contrasts)))
      stop("'contrasts' must be a named vector")

    if (is.null(group.col))
      stop("Must specify a 'group.col' when providing a 'contrasts' argument.")
    
    lmbds = contrasts

    ## Assert that contrasts sum to 0
    if (!isTRUE(all.equal(sum(lmbds), 0, tol=1e-2)))
      stop("'contrasts' must sum to 0.")
    
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
    dat  = dat[dat[[group.col]] %in% boot.groups, ]

    vals = dat[[data.col]]
    grps = as.factor(dat[[group.col]])
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
    if (is.null(group.col))
      stop("You must specify a grouping column.")

    if (length(unique(grps)) != 2)
      stop("There must be only 2 distinct groups!")

    ## Assert that there are 2 numeric columns and reorder these
    ## to be the first 2 columns in the data frame
    num.col.idx = which(sapply(dat, is.numeric))
    num.col.idx = num.col.idx[!names(num.col.idx) %in% group.col]    
    is.valid = length(num.col.idx) == 2
    if (!is.valid)
      stop("'dat' must contain 2 numeric columns and a grouping column.")
    
    group.col.idx = match(group.col, names(dat))
    dat = dat[, c(num.col.idx, group.col.idx)]
  }    
  
  ## Error handling for 'glass.control'
  if (!is.null(glass.control)) {
    if (!is.character(glass.control) || length(glass.control) != 1)
      stop("glass.control must be a character vector of length 1.")

    if (!glass.control %in% colnames(dat))
      stop("'glass.control' missing from 'dat'")
  }  

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
      res = boot(vals, calcSlope, R=R, stype="f", strata=grps,
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
        "with data.col of class %s and group.col of class %s not implemented.",
        collapse="")
      msg = sprintf(msg, effect.type, class(dat), length(dat), class(grps),
        class(vals))
      stop(msg)
    }
  }

  ## Calculate the confidence interval
  if (ci.type != "none") {
    ci = boot.ci(res, conf=ci.conf, type=ci.type, ...)
    ci = ci[[ci.type, exact=FALSE]]

    bounds = switch(ci.type,
      norm  = ci[1, 2:3, drop=TRUE],
      ci[1, 4:5, drop=TRUE])
  } else {
    bounds = c(NA_real_, NA_real_)
  }

  res[["bounds"]]  = bounds
  res[["ci.type"]] = ci.type
  res[["ci.conf"]] = ci.conf
  res[["contrasts"]] = lmbds.orig
  res[["contrasts.scaled"]] = lmbds
  class(res) = c("bootES", "boot")
  
  return(res)
}

print.bootES <- function(x, ...) {
  ## Prints the bootstrap summary statistics for a bootES object, followed by
  ## the confidence interval if the user required its calculation.
  printTerse(x)
  if (isTRUE(x$plot))
    plot(x)
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

  ## BEGIN: Code from boot::print.boot
  index  = seq_len(ncol(x$t))
  t      = matrix(x$t[, index], nrow = nrow(x$t))
  allNA  = apply(t, 2L, function(t) all(is.na(t)))
  ind1   = index[allNA]
  index  = index[!allNA]
  t      = matrix(t[, !allNA], nrow = nrow(t))
  t0     = x$t0

  bias = apply(t, 2L, mean, na.rm=TRUE) - t0
  std.error = sqrt(apply(t, 2L, function(t.st) var(t.st[!is.na(t.st)])))
  ## END: Code from boot::print.boot
  
  nms = c("Stat", "CI (Low)", "CI (High)", "bias", "std. error")
  res = matrix(c(t0, x[["bounds"]], bias, std.error),
               nrow=1, dimnames=list(NULL, nms))

  cat(sprintf("%.2f%% %s Confidence Interval, %d replicates\n",
              100 * x[["ci.conf"]],
              x[["ci.type"]],
              x[["R"]]))
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
  ##  data.col
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
    if (is.null(data.col) && do.cor) {
      res = 'cor'
    } else if (identical(effect.type, 'r')) {
      res = 'r'
    } else {
      res = 'mean'
    }
  } else { 
    if (identical(effect.type, 'slope')) {
      res = 'slope'
    } else if (!is.null(contrasts)) {
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
