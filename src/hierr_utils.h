#ifndef HIERR_UTILS_H
#define HIERR_UTILS_H

#include <RcppArmadillo.h>

Rcpp::NumericVector compute_penalty(Rcpp::NumericVector & ulam,
                                    const int & nlam,
                                    const double & ptype,
                                    const double & pratio,
                                    arma::vec & g,
                                    const arma::vec & cmult,
                                    const int & start,
                                    const int & stop);

void compute_coef(arma::mat & coef,
                  const arma::mat & ext_,
                  const int & nvar,
                  const int & nvar_ext,
                  const int & nvar_total,
                  const int & nlam_total,
                  const arma::vec & xm,
                  const arma::vec & xs,
                  const double & ys,
                  const arma::vec & a0,
                  const bool & intr_ext,
                  const int & ext_start);

void compute_coef_sparse(arma::mat & coef,
                         const arma::sp_mat & ext_,
                         const int & nvar,
                         const int & nvar_ext,
                         const int & nvar_total,
                         const int & nlam_total,
                         const arma::vec & xm,
                         const arma::vec & xs,
                         const double & ys,
                         const arma::vec & a0,
                         const bool & intr_ext,
                         const int & ext_start);

void standardize_vec(arma::vec & y,
                     const arma::vec & w,
                     double & ym,
                     double & ys,
                     const bool & intr);

#endif
