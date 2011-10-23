\name{bootES}
\alias{bootES}
\title{bootES}
\description{
  Provides an interface to various methods for calculating effect sizes using the 'boot' package.
}
\usage{
bootES(dat, R, data.col, grp.col, 
       stat = c("mean", "contrast", "cor", "cor.diff", "slope"), 
       effect.type = c("unstandardized", "cohens.d", "hedges.g", "cohens.d.sigma", "r"),
       contrasts,
       glass.control,
       scale.weights = FALSE,
       ci.type = c("bca", "norm", "basic", "stud", "perc", "all", "none"),
       ci.conf = 0.95)   
}

\arguments{
  \item{dat}{ a numeric vector or data frame containing one or more columns of values
  	      (required), group labels (optional), and contrasts (optional)
              for each sample}
  \item{R}{the number of bootstrap repetitions to perform}
  \item{data.col}{the name of the column in \code{dat} containing the sample values}
  \item{grp.col}{the name of the column in \code{dat} containing the grouping info}
  \item{stat}{The statistic to calculate}
  \item{effect.type}{The unit to report the statistic in}
  \item{contrasts}{An optional vector of contrast weights}
  \item{glass.control}{An optional name of a control group to use}
  \item{scale.weights}{TRUE/FALSE, scale the lambdas to [-1, 1]}
  \item{ci.type}{The method to use in constructing the confidence
    interval.}
  \item{ci.conf}{The confidence level for the confidence interval.}
  \item{verbose}{The amount of output the function generates as it
    runs. May be 0, 1, or 2. When 0, calling 'bootES' does not produce
    any output. When 1, outputs the statistic and confidence
    interval. When 2, outputs all the output from 1 plus diagnostic plots.}
}

\details{Fill this in...}

\value{An object of class \code{bootES}}
\author{Daniel Gerlanc and Kris Kirby}

