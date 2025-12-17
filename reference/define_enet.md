# Define elastic net regularization object for predictor and external data

Helper function to define a elastic net penalty regularization object.
See `define_penalty` for more details.

## Usage

``` r
define_enet(
  en_param = 0.5,
  num_penalty = 20,
  penalty_ratio = NULL,
  user_penalty = NULL,
  custom_multiplier = NULL
)
```

## Arguments

- en_param:

  elastic net parameter, between 0 and 1

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
[`tune_xrnet`](https://uscbiostats.github.io/xrnet/reference/tune_xrnet.md).
The list elements will match those returned by
[`define_penalty`](https://uscbiostats.github.io/xrnet/reference/define_penalty.md),
but with the penalty_type set to match the value of `en_param`.
