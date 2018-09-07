#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * Computes penalty path for set of variables
 * based on simple least-squares estimates, penalty
 * type, number of penalties, and ratio between
 * max and min penalty
 */

//[[Rcpp:export]]
NumericVector compute_penalty(NumericVector & ulam,
                              const int & nlam,
                              const double & ptype,
                              const double & pratio,
                              arma::vec & g,
                              const arma::vec & cmult,
                              const int & start,
                              const int & stop) {

    NumericVector lambdas;
    if (ulam[0] == 0.0) {
        lambdas = NumericVector(nlam);
        lambdas[0] = 9.9e35;
        double max_pen = 0.0;
        for (int j = start; j < stop; ++j) {
            if (cmult[j] > 0.0) {
                max_pen = std::max(max_pen, std::abs(g[j]) / cmult[j]);
            }
        }
        double eqs = std::max(1e-6, pratio);
        double alf = pow(eqs, 1.0 / (nlam - 1));
        lambdas[1] = alf * (max_pen / (std::max(ptype, 0.001)));
        for (int l = 2; l < nlam; l++) {
            lambdas[l] = alf * lambdas[l - 1];
        }
    } else {
        lambdas = ulam;
    }
    return(lambdas);
}

/*
 * Compute and unstandardize estimates
 * for predictor variables (x) as:
 * b = ext * alpha + gamma
 */

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
                  const int & ext_start) {

    if (nvar_ext > 0) {
        for (int j = 0; j < nlam_total; ++j) {
            for (int i = 0; i < nvar; ++i) {
                double z_alpha = a0[j];
                for (int k = ext_start; k < nvar_total; ++k) {
                    z_alpha += coef.at(k, j) * (ext_.at(i, k - ext_start) - xm[k]) / xs[k];
                }
                coef.at(i, j) = ys * (z_alpha + coef.at(i, j)) / xs[i];
            }
        }
    } else {
        for (int j = 0; j < nlam_total; ++j) {
            for (int i = 0; i < nvar; ++i) {
                coef.at(i, j) = ys * coef.at(i, j) / xs[i];
            }
        }
    }
}

/*
 * Compute and unstandardize estimates
 * for predictor variables (x) as:
 * b = ext * alpha + gamma
 */

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
                         const int & ext_start) {

    if (nvar_ext > 0) {
        for (int j = 0; j < nlam_total; ++j) {
            for (int i = 0; i < nvar; ++i) {
                double z_alpha = a0[j];
                for (int k = ext_start; k < nvar_total; ++k) {
                    z_alpha += coef.at(k, j) * (ext_.at(i, k - ext_start) - xm[k]) / xs[k];
                }
                coef.at(i, j) = ys * (z_alpha + coef.at(i, j)) / xs[i];
            }
        }
    } else {
        for (int j = 0; j < nlam_total; ++j) {
            for (int i = 0; i < nvar; ++i) {
                coef.at(i, j) = ys * coef.at(i, j) / xs[i];
            }
        }
    }
}

/*
 * Generic function to standardize vector by
 * mean and standard deviation (using 1/n)
 */

// [[Rcpp::export]]
void standardize_vec(arma::vec & y,
                     const arma::vec & w,
                     double & ym,
                     double & ys,
                     const bool & intr) {
    if (intr) {
        ym = arma::dot(y, w);
        ys = sqrt(arma::dot(y - ym, (y - ym) % w));
        y = (y - ym) / ys;
    } else {
        double ym_temp = arma::dot(w, y);
        ys = sqrt(arma::dot(y - ym_temp, (y - ym_temp) % w));
        y = y / ys;
        ym = 0.0;
    }
}
