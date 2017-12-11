#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}
<<<<<<< HEAD:src/timestwo.cpp
=======

>>>>>>> f53ac40f66e3d39abbaa039b299563f44a44a086:src/timesTwo.cpp
