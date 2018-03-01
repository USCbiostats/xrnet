#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * Coordinate descent function to fit model
 * for design matrix [x | x*ext] on y given
 * variable-specific penalty values and penalty types
 */

// [[Rcpp::export]]
void coord_desc(const arma::mat & x,
                arma::vec & resid,
                const arma::vec & w,
                const NumericVector & ptype,
                const double & tau,
                const double & tau_ext,
                const int & no,
                const int & nvar,
                const int & nvar_total,
                const arma::vec & cmult,
                const NumericVector & upper_cl,
                const NumericVector & lower_cl,
                const int & ne,
                const int & nx,
                NumericVector cur_lam,
                const double & thr,
                const int & maxit,
                const arma::vec & xv,
                arma::mat & coef,
                arma::vec & b,
                arma::vec & g,
                NumericVector & dev,
                double & dev_cur,
                IntegerVector & mm,
                double & errcode,
                int & nlp,
                int & idx_lam) {

    int nin = 0;
    bool jerr = 0;
    NumericVector lasso_part(nvar_total);
    NumericVector lasso_part2(nvar_total);
    NumericVector ridge_part(nvar_total);

    for (int k = 0; k < nvar; k++) {
        lasso_part[k] = 2 * cmult[k] * ptype[k] * cur_lam[0] * tau;
        lasso_part2[k] = 2 * ptype[k] * cur_lam[0];
        ridge_part[k] = xv[k] + cmult[k] * (1 - ptype[k]) * cur_lam[0];
    }

    for (int k = nvar; k < nvar_total; k++) {
        lasso_part[k] = 2 * cmult[k] * ptype[k] * cur_lam[1] * tau_ext;
        lasso_part2[k] = 2 * ptype[k] * cur_lam[1];
        ridge_part[k] = xv[k] + cmult[k] * (1 - ptype[k]) * cur_lam[1];
    }

    LogicalVector active(nvar_total, 1);
    bool last_loop = false;

    bool kkt_satisfied = 0, converge = 0;
    while(!kkt_satisfied & !jerr) {
        while(!converge & !jerr) {
            double dlx = 0.0;
            for (int k = 0; k < nvar_total; ++k) {
                if (active[k]) {
                    double gk = arma::dot(x.unsafe_col(k), resid % w);
                    double bk = b[k];
                    double u = gk + bk * xv[k];
                    double v = std::abs(u) - ((u > 0.0) ? 1 : -1) * lasso_part[k] - lasso_part2[k] * (u < 0.0);
                    if (v > 0.0) {
                        b[k] = std::max(lower_cl[k], std::min(upper_cl[k], copysign(v, u) / ridge_part[k]));
                    } else {
                        b[k] = 0.0;
                    }
                    if (std::abs(b[k] - bk) > 1e-15) {
                        if (mm[k] == 0) {
                            nin += 1;
                            mm[k] = nin;
                        }
                        if (nin <= nx) {
                            double del = b[k] - bk;
                            resid -= del * x.unsafe_col(k);
                            dev_cur += del * (2.0 * gk - del * xv[k]);
                            dlx = std::max(xv[k] * del * del, dlx);
                        } else {
                            jerr = 1;
                            errcode = -10000 - k;
                        }
                    }
                    else {
                        active[k] = 0;
                    }
                }
            }
            nlp += 1;
            if (nlp > maxit) {
                jerr = 1;
                errcode = -10000;
            }
            if (dlx < thr) {
                if (last_loop) {
                    converge = 1;
                }
                else {
                    std::fill(active.begin(), active.end(), 1);
                    last_loop = true;
                }
            }
            else if (last_loop) {
                last_loop = false;
            }
        }
        kkt_satisfied = 1;
    }
    coef.col(idx_lam) = b;
    dev[idx_lam] = dev_cur;
}
