## Submission

This is a major release that fixes gcc-UBSAN errors recently detected in r-devel and changes the minimum R version to be >= 4.0.

* Checks reran on all other platforms as well Rdevel to verify gcc-UBSAN errors were fixed for the case where we try to create a mapped matrix via a pointer to and empty matrix. The mapped matrix is now correctly created in the case the input is an empty matrix.
