#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "hierr_utils.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


/*
 * Update lambdas
 */
void updatePenalty(arma::vec & l1,
                   arma::vec & l2,
                   const arma::vec & pind,
                   const arma::vec & cmult,
                   const arma::vec & xv,
                   const double & lam_new,
                   const int & start,
                   const int & stop) {
    for (int k = start; k < stop; ++k) {
        l1[k] = pind[k] * lam_new;
        l2[k] = 1 / (xv[k] + (cmult[k] - pind[k]) * lam_new);
    }
}

/*
 * Check strong rules
 */

void updateStrong(LogicalVector & strong,
                 const arma::vec & g,
                 const arma::vec & ptype,
                 const arma::vec & cmult,
                 const NumericVector & lam_cur,
                 const NumericVector & lam_prev,
                 const NumericVector & qnt,
                 const IntegerVector & blkend) {
    int nblocks = lam_cur.length();
    int begin = 0;
    for (int blk = 0; blk < nblocks; ++blk) {
        double q = qnt[blk];
        double end = blkend[blk];
        double lam_diff = 2.0 * lam_cur[blk] - lam_prev[blk];
        for (int k = begin; k < end; ++k) {
            if (!strong[k]) {
                strong[k] = std::abs(g[k]) > ptype[k] * cmult[k] * lam_diff * std::abs(q + sgn(g[k]));
            }
        }
        begin = end;
    }
}

/*
 * Computes penalty path for set of variables
 * based on simple least-squares estimates, penalty
 * type, number of penalties, and ratio between
 * max and min penalty
 */

void compute_penalty(NumericVector & lambdas,
                     NumericVector & ulam,
                     const double & ptype,
                     const double & pratio,
                     const arma::vec & g,
                     const arma::vec & cmult,
                     const int & start,
                     const int & stop) {
    int nlam = lambdas.length();
    if (ulam[0] == 0.0) {
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

int countNonzero(const arma::vec & x,
                 const int & start,
                 const int & end) {
    int nzero = 0;
    for(int i = start; i < end; ++i) {
        nzero += (x[i] != 0.0);
    }
    return nzero;
}
