#ifndef COORD_DESC_H
#define COORD_DESC_H

#include <RcppArmadillo.h>

void coord_desc(const arma::mat & x,
                arma::vec & resid,
                const arma::vec & w,
                const arma::vec & ptype_ind,
                const arma::vec & lasso_part,
                const arma::vec & ridge_part,
                const double & qx,
                const double & qext,
                const int & nv_x,
                const int & nvar_total,
                const Rcpp::NumericVector & upper_cl,
                const Rcpp::NumericVector & lower_cl,
                const int & ne,
                const int & nx,
                Rcpp::NumericVector & lam_cur,
                Rcpp::NumericVector & lam_prev,
                Rcpp::LogicalVector & strong,
                Rcpp::IntegerVector & active_x,
                Rcpp::IntegerVector & active_ext,
                const double & thr,
                const int & maxit,
                const arma::vec & xv,
                arma::mat & coef,
                arma::vec & b,
                arma::vec & g,
                Rcpp::NumericVector & dev,
                double & dev_cur,
                Rcpp::LogicalVector & ever_active,
                double & errcode,
                int & nlp,
                int & idx_lam,
                int & nin_x,
                int & nin_ext);
#endif
