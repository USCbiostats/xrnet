# Define regularization object for predictor and external data.

Defines regularization for predictors and external data variables in
[`xrnet`](https://uscbiostats.github.io/xrnet/reference/xrnet.md)
fitting. Use helper functions define_lasso, define_ridge, or define_enet
to specify a common penalty on x or external.

## Usage

``` r
define_penalty(
  penalty_type = 1,
  quantile = 0.5,
  num_penalty = 20,
  penalty_ratio = NULL,
  user_penalty = NULL,
  custom_multiplier = NULL
)
```

## Arguments

- penalty_type:

  type of regularization. Default is 1 (Lasso). Can supply either a
  scalar value or vector with length equal to the number of variables
  the matrix.

  - 0 = Ridge

  - (0,1) = Elastic-Net

  - 1 = Lasso / Quantile

- quantile:

  specifies quantile for quantile penalty. Default of 0.5 reduces to
  lasso (currently not implemented).

- num_penalty:

  number of penalty values to fit in grid. Default is 20.

- penalty_ratio:

  ratio between minimum and maximum penalty for x. Default is 1e-04 if
  \\n \> p\\ and 0.01 if \\n \<= p\\.

- user_penalty:

  user-defined vector of penalty values to use in penalty path.

- custom_multiplier:

  variable-specific penalty multipliers to apply to overall penalty.
  Default is 1 for all variables. 0 is no penalization.

## Value

A list object with regularization settings that are used to define the
regularization for predictors or external data in
[`xrnet`](https://uscbiostats.github.io/xrnet/reference/xrnet.md) and
[`tune_xrnet`](https://uscbiostats.github.io/xrnet/reference/tune_xrnet.md):

- penalty_type:

  The penalty type, scalar with value in range \[0, 1\].

- quantile:

  Quantile for quantile penalty, 0.5 defaults to lasso (not currently
  implemented).

- num_penalty:

  The number of penalty values in the penalty path.

- penalty_ratio:

  The ratio of the minimum penalty value compared to the maximum penalty
  value.

- user_penalty:

  User-defined numeric vector of penalty values, NULL if not provided by
  user.

- custom_multiplier:

  User-defined feature-specific penalty multipliers, NULL if not
  provided by user.

## Examples

``` r
# define ridge penalty with penalty grid split into 30 values
my_penalty <- define_penalty(penalty_type = 0, num_penalty = 30)

# define elastic net (0.5) penalty with user-defined penalty
my_custom_penalty <- define_penalty(
  penalty_type = 0.5, user_penalty = c(100, 50, 10, 1, 0.1)
)
```
