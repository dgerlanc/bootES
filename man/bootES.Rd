\name{bootES}
\alias{bootES}
\title{bootES}
\description{
  Provides an interface to various methods for calculating effect sizes
  using the 'boot' package.
}
\usage{
bootES(dat, R=1000, data.col=NULL, group.col=NULL, 
       effect.type = c("unstandardized", "cohens.d", "hedges.g",
                       "cohens.d.sigma", "r", "slope"),
       contrasts = NULL,
       slope.levels = NULL,
       glass.control = NULL,
       scale.weights = TRUE,
       ci.type = c("bca", "norm", "basic", "perc", "stud", "none"),
       ci.conf = 0.95,
       plot = FALSE,
       ...)   
}

\arguments{
  \item{dat}{a numeric vector or data frame containing one or more columns of values
    (required), group labels (optional), and contrasts (optional)
    for each sample}
  \item{R}{the number of bootstrap repetitions to perform}
  \item{data.col}{the name of the column in \code{dat} containing the sample values}
  \item{group.col}{the name of the column in \code{dat} containing the grouping info}
  \item{effect.type}{The unit to report the statistic in}
  \item{contrasts}{A named vector of between-group contrast weights}
  \item{slope.levels}{A named vector of values corresponding to the
    treatment levels for groups}
  \item{glass.control}{An optional name of a control group to use}
  \item{scale.weights}{TRUE/FALSE, scale the lambdas to [-1, 1]}
  \item{ci.type}{The method to use in constructing the confidence
    interval.}
  \item{ci.conf}{The confidence level for the confidence interval.}
  \item{plot}{Generate the histogram of the sample statistics and the
    quantile-quantile plot on printing?}
  \item{...}{Additional arguments passed to \code{boot.ci}}
}

\value{An object of class \code{bootES}}
\author{Daniel Gerlanc and Kris Kirby}

