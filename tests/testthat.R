Sys.setenv("R_TESTS" = "")
library(testthat)
library(xrnet)
test_check("xrnet")
