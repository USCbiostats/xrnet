# Miscelleneous utility functions used in exported functions

# Function to define loss function for k-fold cross-validation
error_match <- function(family, type.measure) {
    if (family == "gaussian") {
        if(type.measure == "default") {
            type.measure <- "mse"
        }
        else if(!(type.measure %in% c("mse", "deviance", "mae"))) {
            stop(paste("type.measure = ",type.measure," does not match available choices: mse, deviance, mae", sep = ""))
        }
        errfunc <- switch(type.measure,
                      mse = function(beta, y, x, wgt) {wgt * (y - x %*% beta)^2},
                      deviance = function(beta, y, x, wgt) {wgt * (y - x %*% beta)^2},
                      mae = function(beta, y, x, wgt) {wgt * abs(y - x %*% beta)})
    }
    return(errfunc)
}
