<div align="center">

<img src="./.graphics/logo.png" alt="bootES Logo" width=250></img>

# bootES

`bootES` ("bootstrap Effect Sizes") calculates bootstrap confidence intervals for
(un)standardized effect-size measures using the bootstrap.

<!-- badges: start -->
[![CRAN
status](https://www.r-pkg.org/badges/version/bootES)](https://cran.r-project.org/package=bootES)
[![R-CMD-check](https://github.com/dgerlanc/bootes/workflows/R-CMD-check/badge.svg)](https://github.com/dgerlanc/bootes/actions)
<!-- badges: end -->

[Overview](#overview) •
[Install](#install) •
[Get help](#get-help) •
[Contribute](#contribute)

</div>

## Overview

`bootES` uses the 'boot' package to find bootstrap confidence intervals for
unstandardized and standardized effect-size measures appropriate for
experimental and survey research.

Calculate effect sizes for:

- mean effects
- mean differences
- contrasts
- correlations
- differences between correlations

A pre-print version of our article, _bootES:_ An R Package for Bootstrap
Confidence Intervals on Effect Sizes, is available [here](https://web.williams.edu/Psychology/Faculty/Kirby/bootes-kirby-gerlanc-in-press.pdf).

## Install

``` r
# Install from CRAN
install.packages("bootES")

# Install from Github
# install.packages("devtools")
devtools::install_github("dgerlanc/bootES")
```

## Get help
Report bugs by opening an [issue][issues]. If you have a question regarding the
usage of `bootES`, start a [discussion][discussions].

## Contributing

Issues may be filed using [Github Issues](https://github.com/dgerlanc/bootES/issues).

[issues]: https://github.com/dgerlanc/bootES/issues
[discussions]: https://github.com/dgerlanc/bootES/discussions
