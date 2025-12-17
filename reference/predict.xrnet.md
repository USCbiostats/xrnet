# Predict function for "xrnet" object

Extract coefficients or predict response in new data using fitted model
from an
[`xrnet`](https://uscbiostats.github.io/xrnet/reference/xrnet.md)
object. Note that we currently only support returning coefficient
estimates that are in the original path(s).

## Usage

``` r
# S3 method for class 'xrnet'
predict(
  object,
  newdata = NULL,
  newdata_fixed = NULL,
  p = NULL,
  pext = NULL,
  type = c("response", "link", "coefficients"),
  ...
)
```

## Arguments

- object:

  A [`xrnet`](https://uscbiostats.github.io/xrnet/reference/xrnet.md)
  object

- newdata:

  matrix with new values for penalized variables

- newdata_fixed:

  matrix with new values for unpenalized variables

- p:

  vector of penalty values to apply to predictor variables

- pext:

  vector of penalty values to apply to external data variables

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

fit_xrnet <- xrnet(
  x = x_linear,
  y = y_linear,
  external = ext_linear,
  family = "gaussian"
)

lambda1 <- fit_xrnet$penalty[10]
lambda2 <- fit_xrnet$penalty_ext[10]

coef_xrnet <- predict(
  fit_xrnet,
  p = lambda1,
  pext = lambda2,
  type = "coefficients"
)

pred_xrnet <- predict(
  fit_xrnet,
  p = lambda1,
  pext = lambda2,
  newdata = x_linear,
  type = "response"
)
```
