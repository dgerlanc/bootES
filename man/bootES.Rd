\name{bootES}
\alias{bootES}
\title{bootES}
\description{
  \code{bootES} ("bootstrap Effect Sizes") uses the 'boot' package to find bootstrap confidence intervals for unstandardized and standardized effect-size measures appropriate for experimental and survey research. These include effect sizes for mean effects, mean differences, contrasts, correlations, and differences between correlations.
  
}
\usage{
bootES(data, R=2000, data.col=NULL, group.col=NULL, block.col=NULL, 
       effect.type = c("unstandardized", "cohens.d", "hedges.g",
                       "cohens.d.sigma", "r"),
       contrast = NULL,
       slope.levels = NULL,
       glass.control = NULL,
       scale.weights = TRUE,
       ci.type = c("bca", "norm", "basic", "perc", "stud", "none"),
       ci.conf = 0.95,
       plot = FALSE,
       ...)   
}

\arguments{
  \item{data}{A numeric vector, or data frame containing at least one column of numerical values, and (optionally) one or more columns of categorical group labels.}
  \item{R}{The number of bootstrap resamples to perform.}
  \item{data.col}{The name of the column in \code{data} containing the sample values (the "dependent" or "outcome" variable).}
  \item{group.col}{The name of the column in \code{data} containing the categorical group labels (the "independent" or "categorical predictor" variable).}
  \item{block.col}{The name of the column in \code{data} containing the 
  categorical group labels for blocking.}
  \item{effect.type}{Specifies the type of effect-size measure to estimate. The default is "unstandardized".}
  \item{contrast}{A named vector that pairs group labels with their contrast weights. See below for usage.}
  \item{scale.weights}{When TRUE (the default), the contrast weights are scaled so that the unstandardized contrast expresses the effect in the original units of measurement. When FALSE, the magnitude of an unstandardized contrast is a function of the scale of the contrast weights.}
  \item{slope.levels}{A named vector that either contains a column of numerical values, or that pairs group labels in \code{group.col} with their numerical values. This expresses the relationship between two numerical variables as a slope, in the original units of measurement. See below for usage.}
  \item{glass.control}{For Cohen's d-type effect sizes, this option specifies a "control group" from which the standard deviation is estimated.}
  \item{ci.type}{The bootstrap method to use in constructing the confidence interval. The default is "bias-corrected and accelerated" ("bca"). Note that to calculate studentized confidence intervals (\code{ci.type="stud"}), the user must pass \code{var.t} and \code{var.t0} arguments.  See the help page for \code{boot.ci} for more information.}
  \item{ci.conf}{The confidence level for the confidence interval. The default is 0.95.}
  \item{plot}{When TRUE, this generates a histogram and Normal quantile-quantile plot of the resampled statistics. See the help page for \code{boot} for more information.}
  \item{...}{Additional arguments may be passed to \code{boot.ci}.}
}

\details{
	An introduction to \code{bootES}, with illustrations of its usage, can be found in:
	Kirby, K. N., & Gerlanc, D. (2012).  BootES: An R Package for Bootstrap Confidence Intervals on Effect Sizes.  Unpublished manuscript.
	
	When \code{data} is a single column of numbers, and no grouping variable is provided, \code{bootES} finds the mean of that column and its bootstrap confidence interval.
	
	The standardized effect-size options are: 
	"r", which expresses the effect as Pearson's r;
	"cohens.d" or "hedges.g", both of which express an effect in standard deviation units, where the standard deviation is computed as the sample estimate of the population standard deviation, e.g., sqrt(SS/(n-1)); and
	"cohens.d.sigma", which expresses an effect in standard deviation units, where the standard deviation of the sample is taken as the population standard deviation, e.g., sqrt(SS/n).
When one of the Cohen's d-type effect sizes are used, the option \code{glass.control = "controlGroupName"} can be used to specify a single group in \code{group.col} that is to be used for computing the standard deviation.
	
	When both a \code{data.col} and a \code{group.col} is provided, \code{bootES} computes the between-subjects contrast specified in the \code{contrast} option, in the format \code{contrast = c(...)}.  Within the parentheses each group label that is involved in the contrast is followed by an equal sign and its contrast weight, with groups separated by commas.  For example, to find the difference between the means of females and males one could use the option \code{contrast = c(female=1, male=-1)}.  As a convenience, when only two groups are involved in the contrast, the user may omit the contrast weights but add single or double quotes: for example, \code{contrast = c("female", "male")}. So, for a data frame called "myDataFrame" containing a column called "DV" and a column called "Gender," the gender difference in the DV, and its bootstrap confidence interval, can be found by
	
\code{bootES(myDataFrame, data.col="DV", group.col="Gender", contrast=c(female=1, male=-1))}
	
	Within-subjects contrasts must be computed in R as a column of contrast scores, prior to using \code{bootES}. The \code{bootES} function can then be used to find the confidence interval for the mean of that column of contrast scores, or to perform between-subjects contrasts on that column of contrast scores.
	
	The relationship between a numerical outcome variable and a predictor variable (either categorical or continuous) can be expresses as a slope using the \code{slope.levels} option. When a column in the data frame contains the values of the predictor variable, the format is \code{slope.levels = "columnName"}.  For example, for the slope of Income as a function of Age:
	
\code{bootES(myDataFrame, data.col="Income", slope.levels="Age")}

When the predictor column contains group labels (e.g., "child", "teen", and "adult") that correspond to numerical values (e.g., 5, 15, and 25 years old, respectively), the pairing of labels and values can be specified as in the \code{contrast} option.  Note that in this case the \code{group.col} argument is required:

\code{bootES(myDataFrame, data.col="Income", group.col="Age", slope.levels=c(child=5, teen=15, adult=25))}  

	When \code{data} contains only a pair of columns of numbers, \code{bootES} finds the Pearson product-moment correlation between those columns and its bootstrap confidence interval.  For example, in the data frame "myDataFrame," the correlation between Income and Age could be found by
	
\code{bootES(myDataFrame[c("Income","Age")])}

When a third column with two groups is included, and that column is specified as a grouping column, \code{bootES} will find the group difference in the correlation. For example, to find the Gender difference in the Income-Age correlation, we could use

\code{bootES(myDataFrame[c("Income","Age","Gender")], group.col = "Gender")}
	
	\code{bootES} provides functions for calculating several common effect-size measures, and a somewhat simpler interface to the \code{boot} and \code{boot.ci} functions in the \code{boot} package. More complex and custom bootstrap analyses may be done by using the \code{boot} function directly.

  
}

\value{An object of class \code{c('bootES', 'boot')}}
\author{Daniel Gerlanc and Kris Kirby}

