#ifndef CREATE_DATA_H
#define CREATE_DATA_H

#include <RcppArmadillo.h>

arma::mat create_data(const int & nobs,
                      const int & nvar,
                      const int & nvar_ext,
                      const int & nvar_unpen,
                      const int & nvar_total,
                      const arma::mat & x,
                      const arma::mat & ext,
                      const arma::mat & unpen,
                      const arma::vec & w,
                      const bool & isd,
                      const bool & isd_ext,
                      const bool & intr,
                      const bool & intr_ext,
                      arma::vec & xm,
                      arma::vec & xv,
                      arma::vec & xs);

#endif
