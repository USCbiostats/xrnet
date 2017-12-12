#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void standarize_mat(const int & nobs,
                    const int & nvar,
                    NumericMatrix & x,
                    NumericVector & w,
                    const bool & isd,
                    const bool & intr,
                    NumericVector & xm,
                    NumericVector & xv,
                    NumericVector & xs) {

    w = w / sum(w);
    NumericVector v = sqrt(w);

    if (intr) {
        for (int j = 0; j < nvar; j++) {
            xm[j] = std::inner_product(x(_, j).begin(), x(_, j).end(), w.begin(), 0.0);
            x(_, j) = x(_, j) - xm[j];
            xv[j] = std::inner_product(x(_, j).begin(), x(_, j).end(), x(_, j).begin(), 0.0) / nobs;
            if (isd) {
                xs[j] = sqrt(xv[j]);
                x(_, j) = x(_, j) / xs[j];
                xv[j] = 1.0;
            } else {
                xs[j] = 1.0;
            }
        }
    } else {
        for (int j = 0; j < nvar; j++) {
            NumericVector xj = v * x(_, j);
            xv[j] = std::inner_product(xj.begin(), xj.end(), xj.begin(), 0.0);
            if (isd) {
                double xm = std::inner_product(v.begin(), v.end(), xj.begin(), 0.0);
                double xm_sqd = xm*xm;
                double vc = xv[j] - xm_sqd;
                xs[j] = sqrt(vc);
                x(_, j) = x(_, j) / xs[j];
                xv[j] = 1.0 + xm_sqd / vc;
            } else {
                xs[j] = 1.0;
            }
        }
    }
}

// **R

/*
void standardize_vec(NumericVector & y,
                     NumericVector & w,
                     double & ym,
                     double & ys,) {
    w = w / sum(w);

    if (intr) {
        ym = std:inner_product(y.begin(), y.end(), w.begin(), 0.0);
        y = y - ym;
        ys = sqrt(std::inner_product(y.begin(), y.end(), y.begin()))
    }
}
*/
