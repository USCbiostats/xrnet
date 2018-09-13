#ifndef COORD_DESC_H
#define COORD_DESC_H

#include <RcppArmadillo.h>

void coord_desc(const arma::mat & x,
                arma::vec & resid,
                const arma::vec & w,
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
                Rcpp::LogicalVector & strong,
                Rcpp::LogicalVector & ever_active,
                const double & thr,
                const int & maxit,
                const arma::vec & xv,
                arma::vec & b,
                arma::vec & g,
                double & dev_cur,
                double & errcode,
                int & nlp,
                int & nin_x,
                int & nin_ext);
#endif
