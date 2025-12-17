# Fit hierarchical regularized regression model

Fits hierarchical regularized regression model that enables the
incorporation of external data for predictor variables. Both the
predictor variables and external data can be regularized by the most
common penalties (lasso, ridge, elastic net). Solutions are computed
across a two-dimensional grid of penalties (a separate penalty path is
computed for the predictors and external variables). Currently support
regularized linear and logistic regression, future extensions to other
outcomes (i.e. Cox regression) will be implemented in the next major
update.

## Usage

``` r
xrnet(
  x,
  y,
  external = NULL,
  unpen = NULL,
  family = c("gaussian", "binomial"),
  penalty_main = define_penalty(),
  penalty_external = define_penalty(),
  weights = NULL,
  standardize = c(TRUE, TRUE),
  intercept = c(TRUE, FALSE),
  control = list()
)
```

## Arguments

- x:

  predictor design matrix of dimension \\n x p\\, matrix options
  include:

  - matrix

  - big.matrix

  - filebacked.big.matrix

  - sparse matrix (dgCMatrix)

- y:

  outcome vector of length \\n\\

- external:

  (optional) external data design matrix of dimension \\p x q\\, matrix
  options include:

  - matrix

  - sparse matrix (dgCMatrix)

- unpen:

  (optional) unpenalized predictor design matrix, matrix options
  include:

  - matrix

- family:

  error distribution for outcome variable, options include:

  - "gaussian"

  - "binomial"

- penalty_main:

  specifies regularization object for x. See
  [`define_penalty`](https://uscbiostats.github.io/xrnet/reference/define_penalty.md)
  for more details.

- penalty_external:

  specifies regularization object for external. See
  [`define_penalty`](https://uscbiostats.github.io/xrnet/reference/define_penalty.md)
  for more details.

- weights:

  optional vector of observation-specific weights. Default is 1 for all
  observations.

- standardize:

  indicates whether x and/or external should be standardized. Default is
  c(TRUE, TRUE).

- intercept:

  indicates whether an intercept term is included for x and/or external.
  Default is c(TRUE, FALSE).

- control:

  specifies xrnet control object. See
  [`xrnet_control`](https://uscbiostats.github.io/xrnet/reference/xrnet_control.md)
  for more details.

## Value

A list of class `xrnet` with components:

- beta0:

  matrix of first-level intercepts indexed by penalty values

- betas:

  3-dimensional array of first-level penalized coefficients indexed by
  penalty values

- gammas:

  3-dimensional array of first-level non-penalized coefficients indexed
  by penalty values

- alpha0:

  matrix of second-level intercepts indexed by penalty values

- alphas:

  3-dimensional array of second-level external data coefficients indexed
  by penalty values

- penalty:

  vector of first-level penalty values

- penalty_ext:

  vector of second-level penalty values

- family:

  error distribution for outcome variable

- num_passes:

  total number of passes over the data in the coordinate descent
  algorithm

- status:

  error status for xrnet fitting

&nbsp;

- 0 = OK

- 1 = Error/Warning

&nbsp;

- error_msg:

  description of error

## Details

This function extends the coordinate descent algorithm of the R package
`glmnet` to allow the type of regularization (i.e. ridge, lasso) to be
feature-specific. This extension is used to enable fitting hierarchical
regularized regression models, where external information for the
predictors can be included in the `external=` argument. In addition,
elements of the R package `biglasso` are utilized to enable the use of
standard R matrices, memory-mapped matrices from the `bigmemory`
package, or sparse matrices from the `Matrix` package.

## References

Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization
Paths for Generalized Linear Models via Coordinate Descent. Journal of
Statistical Software, 33(1), 1-22. URL
http://www.jstatsoft.org/v33/i01/.

Zeng, Y., and Breheny, P. (2017). The biglasso Package: A Memory- and
Computation-Efficient Solver for Lasso Model Fitting with Big Data in R.
arXiv preprint arXiv:1701.05936. URL https://arxiv.org/abs/1701.05936.

Michael J. Kane, John Emerson, Stephen Weston (2013). Scalable
Strategies for Computing with Massive Data. Journal of Statistical
Software, 55(14), 1-19. URL http://www.jstatsoft.org/v55/i14/.

## Examples

``` r
### hierarchical regularized linear regression ###
data(GaussianExample)

## define penalty for predictors and external variables
## default is ridge for predictors and lasso for external
## see define_penalty() function for more details

penMain <- define_penalty(0, num_penalty = 20)
penExt <- define_penalty(1, num_penalty = 20)

## fit model with defined regularization
fit_xrnet <- xrnet(
  x = x_linear,
  y = y_linear,
  external = ext_linear,
  family = "gaussian",
  penalty_main = penMain,
  penalty_external = penExt
)
```
