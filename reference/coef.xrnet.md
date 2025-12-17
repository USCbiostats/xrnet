# Get coefficient estimates from "xrnet" model object.

Returns coefficients from 'xrnet' model. Note that we currently only
support returning coefficient estimates that are in the original
path(s).

## Usage

``` r
# S3 method for class 'xrnet'
coef(object, p = NULL, pext = NULL, ...)
```

## Arguments

- object:

  A [`xrnet`](https://uscbiostats.github.io/xrnet/reference/xrnet.md)
  object.

- p:

  vector of penalty values to apply to predictor variables.

- pext:

  vector of penalty values to apply to external data variables.

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
data(GaussianExample)

fit_xrnet <- xrnet(
  x = x_linear,
  y = y_linear,
  external = ext_linear,
  family = "gaussian"
)

lambda1 <- fit_xrnet$penalty[10]
lambda2 <- fit_xrnet$penalty_ext[10]

coef_xrnet <- coef(
  fit_xrnet,
  p = lambda1,
  pext = lambda2,
)
```
