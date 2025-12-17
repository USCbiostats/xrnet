# Predict function for "tune_xrnet" object

Extract coefficients or predict response in new data using fitted model
from a
[`tune_xrnet`](https://uscbiostats.github.io/xrnet/reference/tune_xrnet.md)
object. Note that we currently only support returning results that are
in the original path(s).

## Usage

``` r
# S3 method for class 'tune_xrnet'
predict(
  object,
  newdata = NULL,
  newdata_fixed = NULL,
  p = "opt",
  pext = "opt",
  type = c("response", "link", "coefficients"),
  ...
)
```

## Arguments

- object:

  A
  [`tune_xrnet`](https://uscbiostats.github.io/xrnet/reference/tune_xrnet.md)
  object

- newdata:

  matrix with new values for penalized variables

- newdata_fixed:

  matrix with new values for unpenalized variables

- p:

  vector of penalty values to apply to predictor variables. Default is
  optimal value in tune_xrnet object.

- pext:

  vector of penalty values to apply to external data variables. Default
  is optimal value in tune_xrnet object.

- type:

  type of prediction to make using the xrnet model, options include:

  - response

  - link (linear predictor)

  - coefficients

- ...:

  pass other arguments to xrnet function (if needed)

## Value

The object returned is based on the value of type as follows:

- response: An array with the response predictions based on the data for
  each penalty combination

- link: An array with linear predictions based on the data for each
  penalty combination

- coefficients: A list with the coefficient estimates for each penalty
  combination. See
  [`coef.xrnet`](https://uscbiostats.github.io/xrnet/reference/coef.xrnet.md).

## Examples

``` r
data(GaussianExample)

## 5-fold cross validation
cv_xrnet <- tune_xrnet(
  x = x_linear,
  y = y_linear,
  external = ext_linear,
  family = "gaussian",
  control = xrnet_control(tolerance = 1e-6)
)

## Get coefficients and predictions at optimal penalty combination
coef_xrnet <- predict(cv_xrnet, type = "coefficients")
pred_xrnet <- predict(cv_xrnet, newdata = x_linear, type = "response")
```
