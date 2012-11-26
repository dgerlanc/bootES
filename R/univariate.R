## Daniel Gerlanc and Kris Kirby (2010-2012)
## Helper functions for bootstrap analyses using the 'boot' package

## ######
## Univariate bootstrap statistics for means
## ######

## Compute the mean of elements 'i' in vector 'v'
meanBoot <- function(v, i)
  mean(v[i])

## Compute the effect size r for a mean effect.
rMean <- function(x)
  mean(x) / sqrt(mean(x)^2 + sum((x-mean(x))^2) / length(x))

## Compute Cohen's d for a mean effect.
dMean <- function(x)
  mean(x) / sd(x)

## Compute the effect size d for a mean effect for resamples in the boot()
## command.
dMeanBoot <- function(x, i)
  mean(x[i]) / sd(x[i])

## Compute Hedge's g for a mean effect
hMean <- function(x) {
  df   = length(x) - 1
  hadj = gamma(df/2) / ( sqrt(df/2)*gamma((df-1)/2) )
  hadj * mean(x) / sd(x)
}

## Compute the effect size g for a mean effect for resamples in the boot()
## command.
hMeanBoot <- function(x, i) {
  vals = x[i]
  df   = length(vals) - 1
  hadj = gamma(df/2) / ( sqrt(df/2) * gamma((df-1)/2))
  hadj * mean(vals) / sd(vals)
}

## Compute the effect size 'Cohen's sigma d' for a mean effect for
## resamples in the boot() command.
dSigmaMeanBoot <- function(x, i)
  mean(x[i]) / sqrt(sum((x[i] - mean(x[i]))^2) / length(x[i]))

## Compute the effect size r for a mean effect for resamples in the boot()
## command.
rMeanBoot <- function(x, i)
  mean(x[i]) / sqrt(mean(x[i])^2 + sum((x[i] - mean(x[i]))^2) / length(x[i]))





