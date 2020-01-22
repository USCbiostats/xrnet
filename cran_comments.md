## Test environments 
* local Windows 10 x64 install, R 3.6.2
* ubuntu 14.04.5 LTS (travis ci), R 3.6.2 / 4.0.0
* OS X (travis ci), R 3.6.2
* windows i386/x64 (appveyor), R 3.6.2 / 4.0.0
* win-builder, R 3.6.2 / 4.0.0

## R CMD check results
There were no ERRORs or WARNINGs

There was 1 NOTE for ubuntu (3.6.2 / 4.0.0):
* checking installed package size ... NOTE
  installed size is 27.7Mb
  sub-directories of 1Mb or more:
    libs  27.2Mb

  The size appears to be due to the use of Rcpp/RcppEigen.

## BH (boost) dependency
Note that when compiling BH package, the following warning appears multiple times:
  warning: dereferencing type-punned pointer will break strict-aliasing rules [-Wstrict-aliasing]

  This appears to be specific to the Boost header files in the BH package and cannot be altered. 
