# Control function for xrnet fitting

Control function for
[`xrnet`](https://uscbiostats.github.io/xrnet/reference/xrnet.md)
fitting.

## Usage

``` r
xrnet_control(
  tolerance = 1e-08,
  max_iterations = 1e+05,
  dfmax = NULL,
  pmax = NULL,
  lower_limits = NULL,
  upper_limits = NULL
)
```

## Arguments

- tolerance:

  positive convergence criterion. Default is 1e-08.

- max_iterations:

  maximum number of iterations to run coordinate gradient descent across
  all penalties before returning an error. Default is 1e+05.

- dfmax:

  maximum number of variables allowed in model. Default is \\ncol(x) +
  ncol(unpen) + ncol(external) + intercept\[1\] + intercept\[2\]\\.

- pmax:

  maximum number of variables with nonzero coefficient estimate. Default
  is \\min(2 \* dfmax + 20, ncol(x) + ncol(unpen) + ncol(external) +
  intercept\[2\])\\.

- lower_limits:

  vector of lower limits for each coefficient. Default is -Inf for all
  variables.

- upper_limits:

  vector of upper limits for each coefficient. Default is Inf for all
  variables.

## Value

A list object with the following components:

- tolerance:

  The coordinate descent stopping criterion.

- dfmax:

  The maximum number of variables that will be allowed in the model.

- pmax:

  The maximum number of variables with nonzero coefficient estimate.

- lower_limits:

  Feature-specific numeric vector of lower bounds for coefficient
  estimates

- upper_limits:

  Feature-specific numeric vector of upper bounds for coefficient
  estimates
