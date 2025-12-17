# xrnet 1.0.1

* Updates author metadata.

# xrnet 1.0.0

* Changes minimum supported R version 4.0.

* Removes C++11 in SystemRequirements since it is guaranteed that R >= 4.0 that C++11 is minimum supported compiler. This also enables removing CXX_STD=CXX11 in Makevars and Makevars.win.

* Fixes gcc-UBSAN errors detected in latest Rdevel checks. Note that these errors do not impact solutions. They occur when no external data matrix is provided and the underlying C++ code attempts to created a mapped matrix with a pointer to an empty matrix. Since the external data matrix is never used in this case, the resulting standard regularized regression with no external data still produces a correct solution.

# xrnet 0.1.7

* Patched release to fix tests on Solaris OS and removed test dependency on glmnet

# xrnet 0.1.6

* First release to CRAN

* Initial release supports linear and logistic hierarchical regularized regression
