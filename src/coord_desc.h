#ifndef COORD_DESC_H
#define COORD_DESC_H

#include <RcppArmadillo.h>

void coord_desc(const arma::mat & x,
                arma::vec & resid,
                const arma::vec & w,
                const arma::vec & ptype,
                const arma::vec & cmult,
                const double & tau,
                const double & tau_ext,
                const int & nvar,
                const int & nvar_total,
                const Rcpp::NumericVector & upper_cl,
                const Rcpp::NumericVector & lower_cl,
                const int & ne,
                const int & nx,
                Rcpp::NumericVector & lam_cur,
                Rcpp::NumericVector & lam_prev,
                Rcpp::LogicalVector & strong,
                std::vector<int> & active,
                const double & thr,
                const int & maxit,
                const arma::vec & xv,
                arma::mat & coef,
                arma::vec & b,
                arma::vec & g,
                Rcpp::NumericVector & dev,
                double & dev_cur,
                Rcpp::IntegerVector & mm,
                double & errcode,
                int & nlp,
                int & idx_lam);
#endif
