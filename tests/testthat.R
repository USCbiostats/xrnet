Sys.setenv("R_TESTS" = "")
library(testthat)
library(hierr)
test_check("hierr")
