context("family check")

test_that("throw error when matching family not found", {

    x <- matrix(runif(10), nrow = 5)
    y <- 1:5
    external <- matrix(runif(10), nrow = 2)
    expect_error(hierr(x, y, external, "badfamily"))
})
