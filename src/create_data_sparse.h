#ifndef CREATE_DATA_SPARSE_H
#define CREATE_DATA_SPARSE_H

#include <RcppArmadillo.h>

arma::mat create_data_sparse(const int & nobs,
                             const int & nvar,
                             const int & nvar_ext,
                             const int & nvar_unpen,
                             const int & nvar_total,
                             const arma::mat & x,
                             const arma::sp_mat & ext,
                             const arma::mat & unpen,
                             const arma::vec & w,
                             const bool & isd,
                             const bool & isd_ext,
                             const bool & intr,
                             const bool & intr_ext,
                             arma::vec & xm,
                             arma::vec & xv,
                             arma::vec & xs,
                             int & ext_start);

double mean_sparse(const arma::sp_mat & x,
                   int & col_cur,
                   const int & nvar);

double sd_sparse(const arma::sp_mat & x,
                 int & col_cur,
                 double & xm,
                 const int & nvar);

#endif
