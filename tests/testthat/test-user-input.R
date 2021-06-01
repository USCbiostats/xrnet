library(bigmemory)
library(Matrix)

context("check user input errors")

######################### xrnet() errors #########################

test_that("throw error when matching family not found", {
  x <- matrix(runif(10), nrow = 5)
  y <- 1:5
  external <- matrix(runif(10), nrow = 2)
  expect_error(xrnet(x, y, external, family = "badfamily"))
})

test_that("throw error when dimensions of x and y do not match", {
  x <- matrix(runif(10), ncol = 2)
  y <- 1:10
  external <- matrix(runif(10), nrow = 2)
  expect_error(
    xrnet(x, y, external, "gaussian"),
    "Length of y (10) not equal to the number of rows of x (5)",
    fixed = TRUE
  )
})

test_that("throw error when dimensions of unpen and y do not match", {
  x <- matrix(runif(20), ncol = 2)
  unpen <- matrix(runif(10), ncol = 2)
  y <- 1:10
  external <- matrix(runif(10), nrow = 2)
  expect_error(
    xrnet(x, y, external, unpen = unpen, family = "gaussian"),
    "Length of y (10) not equal to the number of rows of unpen (5)",
    fixed = TRUE
  )
})

test_that("throw error when ncol(x) not equal to nrow(external)", {
  x <- matrix(runif(10), ncol = 2)
  y <- 1:5
  external <- matrix(runif(10), nrow = 5)
  expect_error(
    xrnet(x, y, external, "gaussian"),
    "Number of columns in x (2) not equal to the number of rows in external (5)",
    fixed = TRUE
  )
})

test_that("throw error if x is not of type double", {
  x <- matrix(1L:10L, nrow = 5)
  y <- 1:5
  expect_error(xrnet(x, y, family = "gaussian"))
  expect_error(xrnet(as.big.matrix(x), y, family = "gaussian"))
})

test_that("throw error if x not one of accepted types", {
  x <- list(1:10)
  y <- 1:5
  expect_error(xrnet(x, y, family = "gaussian"))
  expect_error(xrnet(x, y, family = "gaussian"))
})

test_that("throw error if external is not of type double", {
  x <- matrix(runif(10), nrow = 5)
  y <- 1:5
  external <- matrix(1L:10L, nrow = 2)
  expect_error(xrnet(x, y, external, family = "gaussian"))
})

test_that("throw error when dimensions of weights and y do not match", {
  x <- matrix(runif(10), ncol = 2)
  y <- 1:5
  wgts <- rep(1:3)
  expect_error(
    xrnet(x, y, weights = wgts, family = "gaussian"),
    "Length of weights (3) not equal to length of y (5)",
    fixed = TRUE
  )
})

test_that("throw error when weights negative", {
  x <- matrix(runif(10), ncol = 2)
  y <- 1:5
  wgts <- c(rep(1, 4), -1)
  expect_error(
    xrnet(x, y, weights = wgts, family = "gaussian"),
    fixed = TRUE
  )
})



######################### xrnet_control() errors #########################

test_that("throw error when tolerance non-positive", {
  expect_error(xrnet_control(tolerance = 0))
  expect_error(xrnet_control(tolerance = -1))
})

test_that("throw error when max iterations non-positive or not an integer", {
  expect_error(xrnet_control(max_iterations = 0))
  expect_error(xrnet_control(max_iterations = -1))
  expect_error(xrnet_control(max_iterations = 2.5))
})

######################### initialize_penalty() errors #########################

test_that("throw error when length of penalty_type != ncol(x)", {
  x <- matrix(runif(20), ncol = 5)
  y <- 1:4
  p <- define_penalty(penalty_type = rep(1, 4))
  expect_error(xrnet(x, y, family = "gaussian", penalty_main = p))
})

test_that("throw error when num_penalty < 3", {
  x <- matrix(runif(20), ncol = 5)
  y <- 1:4
  p <- define_penalty(num_penalty = 2)
  expect_error(xrnet(x, y, family = "gaussian", penalty_main = p))
})

test_that("throw error when length of custom_multiplier != ncol(x)", {
  x <- matrix(runif(20), ncol = 5)
  y <- 1:4
  p <- define_penalty(custom_multiplier = rep(10, 4))
  expect_error(xrnet(x, y, family = "gaussian", penalty_main = p))
})

test_that("throw error when length of penalty_type_ext != ncol(external)", {
  x <- matrix(runif(20), ncol = 5)
  y <- 1:4
  external <- matrix(runif(20), nrow = 5)
  p <- define_penalty(penalty_type = rep(0, 2))
  expect_error(xrnet(x, y, external, family = "gaussian", penalty_external = p))
})

test_that("throw error when num_penalty_ext < 3", {
  x <- matrix(runif(20), ncol = 5)
  y <- 1:4
  external <- matrix(runif(20), nrow = 5)
  p <- define_penalty(num_penalty = 2)
  expect_error(xrnet(x, y, external, family = "gaussian", penalty_external = p))
})

test_that("throw error when length of custom_multiplier_ext != ncol(external)", {
  x <- matrix(runif(20), ncol = 5)
  y <- 1:4
  external <- matrix(runif(20), nrow = 5)
  p <- define_penalty(custom_multiplier = rep(10, 2))
  expect_error(xrnet(x, y, external, family = "gaussian", penalty_external = p))
})

######################### initialize_control() errors #########################

test_that("throw error when dfmax non-positive or not an integer", {
  x <- matrix(runif(10), nrow = 5)
  y <- 1:5

  expect_error(
    xrnet(
      x = x,
      y = y,
      family = "gaussian",
      control = xrnet_control(dfmax = 0)
    )
  )

  expect_error(
    xrnet(
      x = x,
      y = y,
      family = "gaussian",
      control = xrnet_control(dfmax = -1)
    )
  )

  expect_error(
    xrnet(
      x = x,
      y = y,
      family = "gaussian",
      control = xrnet_control(dfmax = 2.5)
    )
  )
})


test_that("throw error when pmax non-positive or not an integer", {
  x <- matrix(runif(10), nrow = 5)
  y <- 1:5

  expect_error(
    xrnet(
      x = x,
      y = y,
      family = "gaussian",
      control = xrnet_control(pmax = 0)
    )
  )

  expect_error(
    xrnet(
      x = x,
      y = y,
      family = "gaussian",
      control = xrnet_control(pmax = -1)
    )
  )

  expect_error(
    xrnet(
      x = x,
      y = y,
      family = "gaussian",
      control = xrnet_control(pmax = 2.5)
    )
  )
})

test_that("throw error when length lower_limits or upper_limits does not match total number of variables", {
  x <- matrix(runif(20), ncol = 5)
  y <- 1:5

  expect_error(
    xrnet(
      x = x,
      y = y,
      family = "gaussian",
      control = xrnet_control(lower_limits = rep(0, 2))
    )
  )

  expect_error(
    xrnet(
      x = x,
      y = y,
      family = "gaussian",
      control = xrnet_control(upper_limits = rep(0, 2))
    )
  )
})
