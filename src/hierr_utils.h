#ifndef HIERR_UTILS_H
#define HIERR_UTILS_H

#include <RcppArmadillo.h>

template <class T>
inline T sgn(T & v) {
    return (v > T(0)) - (v < T(0));
}

int countNonzero(const arma::vec & x,
                 const int & start,
                 const int & end);

void updatePenalty(arma::vec & l1,
                   arma::vec & l2,
                   const arma::vec & pind,
                   const arma::vec & cmult,
                   const arma::vec & xv,
                   const double & lam_new,
                   const int & start,
                   const int & stop);

void updateStrong(Rcpp::LogicalVector & strong,
                  const arma::vec & g,
                  const arma::vec & ptype,
                  const arma::vec & cmult,
                  const Rcpp::NumericVector & lam_cur,
                  const Rcpp::NumericVector & lam_prev,
                  const Rcpp::NumericVector & qnt,
                  const Rcpp::IntegerVector & blkend);

void compute_penalty(Rcpp::NumericVector & lambdas,
                     Rcpp::NumericVector & ulam,
                     const double & ptype,
                     const double & pratio,
                     const arma::vec & g,
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
