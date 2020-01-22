## Test environments 
* local Windows 10 x64 install, R 3.6.2
* ubuntu 14.04.5 LTS (travis ci), R 3.6.2 / 4.0.0
* osx (travis ci), R 3.6.2
* windows i386/x64 (appveyor), R 3.6.2 / 4.0.0

## R CMD check results
There were no ERRORs or WARNINGs

There was 1 NOTE for ubuntu (3.6.2 / 4.0.0):
* checking installed package size ... NOTE
  installed size is 27.7Mb
  sub-directories of 1Mb or more:
    libs  27.2Mb

  The size appears to be due to the use of Rcpp/RcppEigen.
