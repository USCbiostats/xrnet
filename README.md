# hierr: An R Package for Hierarchical Regularized Regression

[![Build Status](https://travis-ci.org/gmweaver/hierr.svg?branch=master)](https://travis-ci.org/gmweaver/hierr)
[![codecov](https://codecov.io/gh/gmweaver/hierr/branch/master/graph/badge.svg)](https://codecov.io/gh/gmweaver/hierr)

# Setup

1. Install [RTools](https://cran.r-project.org/bin/windows/Rtools/) (not an R package) and the R package devtools
2. Install hierr 
3. Load the package

```R
library(devtools)
install_github(gmweaver/hierr)
library(hierr)
```

# Getting Started

As an example of how you might use the `hierr` package, we have provided a small simulated set of external data, predictors, and outcome. First, load the example data:

```
load("x.rda")
load("y.rda")
load("ext.rda")
```

To fit a linear hierarchical regularized regression model, you must specify the predictor matrix, outcome, external data, outcome distribution. By default, a ridge penalty is applied to the predictors and a lasso penalty is applied to the external data.

```
hierr_model <- hierr(x = x, y = y, external = ext, family = "gaussian")
```
