# hierr: An R Package for Hierarchical Regularized Regression

[![Build Status](https://travis-ci.org/gmweaver/hierr.svg?branch=master)](https://travis-ci.org/gmweaver/hierr)
[![codecov](https://codecov.io/gh/gmweaver/hierr/branch/master/graph/badge.svg)](https://codecov.io/gh/gmweaver/hierr)

# Setup

1. Install [RTools](https://cran.r-project.org/bin/windows/Rtools/) (not an R package) and the R package devtools
2. Install hierr by using the devtools function install_github
3. Load the package

```R
library(devtools)
install_github(gmweaver/hierr)
library(hierr)
```

# Getting Started

As an example of how you might use the hierr package, we have provided a small example with simulated external data variables (ext), predictors (x), and a continuous outcome variable (y). First, load the example data:

```R
load("x.rda")
load("y.rda")
load("ext.rda")
```

To fit a linear hierarchical regularized regression model, you can use the `hierr` function. At a minimum, you must specify the predictor matrix, outcome variable, external data matrix, and outcome distribution. By default, a ridge penalty is applied to the predictors and a lasso penalty is applied to the external data.

```R
hierr_model <- hierr(x = x, y = y, external = ext, family = "gaussian")
```

To modify the regularization terms and penalty path associated with the predictors or external data, you can use the `definePenalty` function. This function allows you to configure the following regularization attributes for both the predictors (x) and the external data variables (external):

* Regularization type 
    - Ridge = 0
    - Elastic Net = (0, 1)
    - Lasso = 1
* Penalty path
    - Number of penalty values in the full penalty path (default = 20)
    - Ratio of min(penalty) / max(penalty) 
* User-defined set of penalties

As an example, we may want to apply a ridge penalty to both x and external. In addition, we may want to have 30 penalty values computed for the regularization terms associated with both x and external. We modify our model fitting as follows.

```R
myPenalty <- definePenalty(penalty_type = 0, penalty_type_ext = 0, num_penalty = 30, num_penalty_ext = 30)
hierr_model2 <- hierr(x = x, y = y, external = ext, family = "gaussian", penalty = myPenalty)
```
