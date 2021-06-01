#' k-fold cross-validation for hierarchical regularized regression
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom bigmemory describe
#' @importFrom bigmemory attach.big.matrix
#'
#' @description k-fold cross-validation for hierarchical regularized
#' regression \code{\link{xrnet}}
#'
#' @param x predictor design matrix of dimension \eqn{n x p}, matrix options
#' include:
#' \itemize{
#'    \item matrix
#'    \item big.matrix
#'    \item filebacked.big.matrix
#'    \item sparse matrix (dgCMatrix)
#' }
#' @param y outcome vector of length \eqn{n}
#' @param external (optional) external data design matrix of dimension
#' \eqn{p x q}, matrix options include:
#' \itemize{
#'     \item matrix
#'     \item sparse matrix (dgCMatrix)
#' }
#' @param unpen (optional) unpenalized predictor design matrix, matrix options
#' include:
#' \itemize{
#'     \item matrix
#' }
#' @param family error distribution for outcome variable, options include:
#' \itemize{
#'     \item "gaussian"
#'     \item "binomial"
#' }
#' @param penalty_main specifies regularization object for x. See
#' \code{\link{define_penalty}} for more details.
#' @param penalty_external specifies regularization object for external. See
#' \code{\link{define_penalty}} for more details.
#' See \code{\link{define_penalty}} for more details.
#' @param weights optional vector of observation-specific weights.
#' Default is 1 for all observations.
#' @param standardize indicates whether x and/or external should be
#' standardized. Default is c(TRUE, TRUE).
#' @param intercept indicates whether an intercept term is included for x and/or
#' external. Default is c(TRUE, FALSE).
#' @param loss loss function for cross-validation. Options include:
#' \itemize{
#'    \item "deviance"
#'    \item "mse" (Mean Squared Error)
#'    \item "mae" (Mean Absolute Error)
#'    \item "auc" (Area under the curve)
#' }
#' @param nfolds number of folds for cross-validation. Default is 5.
#' @param foldid (optional) vector that identifies user-specified fold for each
#' observation. If NULL, folds are automatically generated.
#' @param parallel use \code{foreach} function to fit folds in parallel if TRUE,
#' must register cluster (\code{doParallel}) before using.
#' @param control specifies xrnet control object. See
#' \code{\link{xrnet_control}} for more details.
#'
#' @return A list of class \code{tune_xrnet} with components
#' \item{cv_mean}{mean cross-validated error for each penalty combination.
#' Object returned is a vector if there is no external data (external = NULL)
#' and matrix if there is external data.}
#' \item{cv_sd}{estimated standard deviation for cross-validated errors.
#' Object returned is a vector if there is no external data (external = NULL)
#' and matrix if there is external data.}
#' \item{loss}{loss function used to compute cross-validation error}
#' \item{opt_loss}{the value of the loss function for the optimal
#' cross-validated error}
#' \item{opt_penalty}{first-level penalty value that achieves the optimal loss}
#' \item{opt_penalty_ext}{second-level penalty value that achieves the optimal
#' loss (if external data is present)}
#' \item{fitted_model}{fitted xrnet object using all data, see
#' \code{\link{xrnet}} for details of object}
#'
#' @details k-fold cross-validation is used to determine the 'optimal'
#' combination of hyperparameter values, where optimal is based on the optimal
#' value obtained for the user-selected loss function across the k folds. To
#' efficiently traverse all possible combinations of the hyperparameter values,
#' 'warm-starts' are used to traverse the penalty from largest to smallest
#' penalty value(s). Note that the penalty grid for the folds is generated
#' by fitting the model on the entire training data. Parallelization is enabled
#' through the \code{foreach} and \code{doParallel} R packages. To use
#' parallelization, \code{parallel = TRUE}, you must first create the cluster
#' \code{makeCluster} and then register the cluster \code{registerDoParallel}.
#' See the \code{parallel}, \code{foreach}, and/or \code{doParallel} R packages
#' for more details on how to setup parallelization.
#'
#' @examples
#' ## cross validation of hierarchical linear regression model
#' data(GaussianExample)
#'
#' ## 5-fold cross validation
#' cv_xrnet <- tune_xrnet(
#'   x = x_linear,
#'   y = y_linear,
#'   external = ext_linear,
#'   family = "gaussian",
#'   control = xrnet_control(tolerance = 1e-6)
#' )
#'
#' ## contour plot of cross-validated error
#' plot(cv_xrnet)
#' @export
tune_xrnet <- function(x,
                       y,
                       external = NULL,
                       unpen = NULL,
                       family = c("gaussian", "binomial"),
                       penalty_main = define_penalty(),
                       penalty_external = define_penalty(),
                       weights = NULL,
                       standardize = c(TRUE, TRUE),
                       intercept = c(TRUE, FALSE),
                       loss = c("deviance", "mse", "mae", "auc"),
                       nfolds = 5,
                       foldid = NULL,
                       parallel = FALSE,
                       control = list()) {
  # function call
  this_call <- match.call()

  # Check family argument
  family <- match.arg(family)

  # Set measure used to assess model prediction performance
  if (missing(loss)) {
    if (family == "gaussian") {
      loss <- "mse"
    } else if (family == "binomial") {
      loss <- "auc"
    }
  } else {
    loss <- match.arg(loss)
    loss_available <- TRUE
    if (family == "gaussian" && !(loss %in% c("deviance", "mse", "mae"))) {
      loss_available <- FALSE
    } else if (family == "binomial" && !(loss %in% c("deviance", "auc"))) {
      loss_available <- FALSE
    }
    if (!loss_available) {
      stop(
        paste0(
          "loss = '",
          loss,
          "' is not available for family = '",
          family,
          "'"
        )
      )
    }
  }

  # check type of x matrix
  if (is(x, "matrix")) {
    if (!(typeof(x) %in% c("integer", "double"))) {
      stop("x contains non-numeric values")
    }
    mattype_x <- 1
  } else if (is.big.matrix(x)) {
    if (
      !(bigmemory::describe(x)@description$type %in% c("integer", "double"))
    ) {
      stop("x contains non-numeric values")
    }
    mattype_x <- 2
  } else if ("dgCMatrix" %in% class(x)) {
    if (!(typeof(x@x) %in% c("integer", "double"))) {
      stop("x contains non-numeric values")
    }
    mattype_x <- 3
  } else {
    stop(
      "x must be a standard R matrix,
      big.matrix, filebacked.big.matrix, or dgCMatrix"
    )
  }

  # check external type
  is_sparse_ext <- is(external, "sparseMatrix")

  # check y type
  y <- drop(as.numeric(y))

  # Get arguments to tune_xrnet() function and filter for calls to fitting
  xrnet_call <- match.call(expand.dots = TRUE)

  cv_args <- match(
    c("loss", "nfolds", "foldid", "parallel"), names(xrnet_call),
    FALSE
  )

  if (any(cv_args)) {
    xrnet_call <- xrnet_call[-cv_args]
  }
  xrnet_call[[1]] <- as.name("xrnet")

  # Set sample size / weights
  n <- length(y)
  if (is.null(weights)) {
    weights <- rep(1, n)
  }

  # Fit model on all training data
  xrnet_object <- xrnet(
    x = x,
    y = y,
    external = external,
    unpen = unpen,
    family = family,
    weights = weights,
    standardize = standardize,
    intercept = intercept,
    penalty_main = penalty_main,
    penalty_external = penalty_external,
    control = control
  )
  xrnet_object$call <- xrnet_call

  # Check whether fixed and external are empty
  if (is.null(unpen)) {
    unpen <- matrix(vector("numeric", 0), 0, 0)
    nc_unpen <- as.integer(0)
  } else {
    nc_unpen <- NCOL(unpen)
  }
  if (is.null(external)) {
    external <- matrix(vector("numeric", 0), 0, 0)
    nc_ext <- as.integer(0)
  } else {
    nc_ext <- NCOL(external)
  }

  # Prepare penalty and control object for folds
  penalty_main_fold <- penalty_main
  penalty_external_fold <- penalty_external

  penalty_main_fold$user_penalty <- xrnet_object$penalty
  if (is.null(xrnet_object$penalty_ext)) {
    penalty_external_fold$user_penalty <- as.double(0.0)
  } else {
    penalty_external_fold$user_penalty <- xrnet_object$penalty_ext
  }

  penalty_fold <- initialize_penalty(
    penalty_main = penalty_main_fold,
    penalty_external = penalty_external_fold,
    nr_x = NROW(x),
    nc_x = NCOL(x),
    nc_unpen = nc_unpen,
    nr_ext = NROW(external),
    nc_ext = nc_ext,
    intercept = intercept
  )

  num_pen <- penalty_fold$num_penalty
  num_pen_ext <- penalty_fold$num_penalty_ext

  control <- do.call("xrnet_control", control)
  control <- initialize_control(
    control_obj = control,
    nc_x = NCOL(x),
    nc_unpen = nc_unpen,
    nc_ext = nc_ext,
    intercept = intercept
  )

  # Randomly sample observations into folds / check nfolds
  if (is.null(foldid)) {
    if (nfolds < 2) {
      stop("number of folds (nfolds) must be at least 2")
    }
    foldid <- sample(rep(seq(nfolds), length = n))
  } else {
    if (length(foldid) != n) {
      stop(
        "length of foldid (", length(foldid), ")
        not equal to number of observations (", n, ")"
      )
    }
    foldid <- as.numeric(factor(foldid))
    nfolds <- length(unique(foldid))
    if (nfolds < 2) {
      stop("number of folds (nfolds) must be at least 2")
    }
  }

  # Run k-fold CV
  if (parallel) {
    if (is.big.matrix(x)) {
      xdesc <- describe(x)
      errormat <- foreach(
        k = 1L:nfolds,
        .packages = c("bigmemory", "xrnet"),
        .combine = cbind
      ) %dopar% {
        weights_train <- weights
        weights_train[foldid == k] <- 0.0
        test_idx <- as.integer(which(foldid == k) - 1)
        xref <- attach.big.matrix(xdesc)

        error_vec <- fitModelCVRcpp(
          x = xref,
          mattype_x = mattype_x,
          y = y,
          ext = external,
          is_sparse_ext = is_sparse_ext,
          fixed = unpen,
          weights_user = weights_train,
          intr = intercept,
          stnd = standardize,
          penalty_type = penalty_fold$ptype,
          cmult = penalty_fold$cmult,
          quantiles = c(
            penalty_fold$quantile, penalty_fold$quantile_ext
          ),
          num_penalty = c(
            penalty_fold$num_penalty, penalty_fold$num_penalty_ext
          ),
          penalty_ratio = c(
            penalty_fold$penalty_ratio, penalty_fold$penalty_ratio_ext
          ),
          penalty_user = penalty_fold$user_penalty,
          penalty_user_ext = penalty_fold$user_penalty_ext,
          lower_cl = control$lower_limits,
          upper_cl = control$upper_limits,
          family = family,
          user_loss = loss,
          test_idx = test_idx,
          thresh = control$tolerance,
          maxit = control$max_iterations,
          ne = control$dfmax,
          nx = control$pmax
        )
      }
    } else {
      errormat <- foreach(
        k = 1L:nfolds,
        .packages = c("Matrix", "xrnet"),
        .combine = cbind
      ) %dopar% {
        weights_train <- weights
        weights_train[foldid == k] <- 0.0
        test_idx <- as.integer(which(foldid == k) - 1)

        # Get errors for k-th fold
        error_vec <- fitModelCVRcpp(
          x = x,
          mattype_x = mattype_x,
          y = y,
          ext = external,
          is_sparse_ext = is_sparse_ext,
          fixed = unpen,
          weights_user = weights_train,
          intr = intercept,
          stnd = standardize,
          penalty_type = penalty_fold$ptype,
          cmult = penalty_fold$cmult,
          quantiles = c(
            penalty_fold$quantile, penalty_fold$quantile_ext
          ),
          num_penalty = c(
            penalty_fold$num_penalty, penalty_fold$num_penalty_ext
          ),
          penalty_ratio = c(
            penalty_fold$penalty_ratio, penalty_fold$penalty_ratio_ext
          ),
          penalty_user = penalty_fold$user_penalty,
          penalty_user_ext = penalty_fold$user_penalty_ext,
          lower_cl = control$lower_limits,
          upper_cl = control$upper_limits,
          family = family,
          user_loss = loss,
          test_idx = test_idx,
          thresh = control$tolerance,
          maxit = control$max_iterations,
          ne = control$dfmax,
          nx = control$pmax
        )
      }
    }
  } else {
    errormat <- matrix(NA, nrow = num_pen * num_pen_ext, ncol = nfolds)
    for (k in 1:nfolds) {
      # Split into test and train for k-th fold
      weights_train <- weights
      weights_train[foldid == k] <- 0.0
      test_idx <- as.integer(which(foldid == k) - 1)

      # Fit model on k-th training fold
      errormat[, k] <- fitModelCVRcpp(
        x = x,
        mattype_x = mattype_x,
        y = y,
        ext = external,
        is_sparse_ext = is_sparse_ext,
        fixed = unpen,
        weights_user = weights_train,
        intr = intercept,
        stnd = standardize,
        penalty_type = penalty_fold$ptype,
        cmult = penalty_fold$cmult,
        quantiles = c(
          penalty_fold$quantile, penalty_fold$quantile_ext
        ),
        num_penalty = c(
          penalty_fold$num_penalty, penalty_fold$num_penalty_ext
        ),
        penalty_ratio = c(
          penalty_fold$penalty_ratio, penalty_fold$penalty_ratio_ext
        ),
        penalty_user = penalty_fold$user_penalty,
        penalty_user_ext = penalty_fold$user_penalty_ext,
        lower_cl = control$lower_limits,
        upper_cl = control$upper_limits,
        family = family,
        user_loss = loss,
        test_idx = test_idx,
        thresh = control$tolerance,
        maxit = control$max_iterations,
        ne = control$dfmax,
        nx = control$pmax
      )
    }
  }
  cv_mean <- rowMeans(errormat)
  cv_sd <- sqrt(rowSums((errormat - cv_mean)^2) / nfolds)
  cv_mean <- matrix(cv_mean, nrow = num_pen, byrow = TRUE)
  cv_sd <- matrix(cv_sd, nrow = num_pen, byrow = TRUE)
  rownames(cv_mean) <- rev(sort(xrnet_object$penalty))
  rownames(cv_sd) <- rev(sort(xrnet_object$penalty))
  if (num_pen_ext > 1) {
    colnames(cv_mean) <- rev(sort(xrnet_object$penalty_ext))
    colnames(cv_sd) <- rev(sort(xrnet_object$penalty_ext))
  }
  if (loss %in% c("deviance", "mse", "mae")) {
    opt_loss <- min(cv_mean, na.rm = TRUE)
    opt_index <- which(opt_loss == cv_mean, arr.ind = TRUE)
  } else {
    opt_loss <- max(cv_mean, na.rm = TRUE)
    opt_index <- which(opt_loss == cv_mean, arr.ind = TRUE)
  }

  if (is.null(dim(opt_index))) {
    opt_penalty <- xrnet_object$penalty[opt_index[1]]
    opt_penalty_ext <- xrnet_object$penalty_ext[opt_index[2]]
  } else {
    opt_penalty <- xrnet_object$penalty[opt_index[1, 1]]
    opt_penalty_ext <- xrnet_object$penalty_ext[opt_index[1, 2]]
  }

  cvfit <- list(
    cv_mean = cv_mean,
    cv_sd = cv_sd,
    loss = loss,
    opt_loss = opt_loss,
    opt_penalty = opt_penalty,
    opt_penalty_ext = opt_penalty_ext,
    fitted_model = xrnet_object,
    call = this_call
  )

  class(cvfit) <- "tune_xrnet"
  return(cvfit)
}
