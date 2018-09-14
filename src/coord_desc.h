#ifndef COORD_DESC_H
#define COORD_DESC_H

#include <RcppArmadillo.h>

void coord_desc(const arma::mat & x,
                arma::vec & resid,
                const arma::vec & w,
                const arma::vec & ptype,
                const arma::vec & cmult,
                const Rcpp::NumericVector & qnt,
                const Rcpp::NumericVector & lam_cur,
                const Rcpp::IntegerVector & stop,
                const Rcpp::NumericVector & ucl,
                const Rcpp::NumericVector & lcl,
                const int & ne,
                const int & nx,
                Rcpp::LogicalVector & strong,
                Rcpp::LogicalVector & active,
                const double & thr,
                const int & maxit,
                const arma::vec & xv,
                arma::vec & b,
                arma::vec & g,
                double & dev_cur,
                double & errcode,
                int & nlp,
                Rcpp::IntegerVector & nin);
#endif
