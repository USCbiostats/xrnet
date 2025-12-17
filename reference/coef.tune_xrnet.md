# Get coefficient estimates from "tune_xrnet" model object.

Returns coefficients from 'xrnet' model. Note that we currently only
support returning coefficient estimates that are in the original
path(s).

## Usage

``` r
# S3 method for class 'tune_xrnet'
coef(object, p = "opt", pext = "opt", ...)
```

## Arguments

- object:

  A
  [`tune_xrnet`](https://uscbiostats.github.io/xrnet/reference/tune_xrnet.md)
  object.

- p:

  vector of penalty values to apply to predictor variables. Default is
  optimal value in tune_xrnet object.

- pext:

  vector of penalty values to apply to external data variables. Default
  is optimal value in tune_xrnet object.

- ...:

  pass other arguments to xrnet function (if needed).

## Value

A list with coefficient estimates at each of the requested penalty
combinations.

- beta0:

  matrix of first-level intercepts indexed by penalty values, NULL if no
  first-level intercept in original model fit.

- betas:

  3-dimensional array of first-level penalized coefficients indexed by
  penalty values.

- gammas:

  3-dimensional array of first-level non-penalized coefficients indexed
  by penalty values, NULL if unpen NULL in original model fit.

- alpha0:

  matrix of second-level intercepts indexed by penalty values, NULL if
  no second-level intercept in original model fit.

- alphas:

  3-dimensional array of second-level external data coefficients indexed
  by penalty values, NULL if external NULL in original model fit.

## Examples

``` r
## Cross validation of hierarchical linear regression model
data(GaussianExample)

## 5-fold cross validation
cv_xrnet <- tune_xrnet(
  x = x_linear,
  y = y_linear,
  external = ext_linear,
  family = "gaussian",
  control = xrnet_control(tolerance = 1e-6)
)

## Get coefficient estimates at optimal penalty combination
coef_opt <- coef(cv_xrnet)
```
