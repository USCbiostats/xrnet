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
                const arma::vec & ptype,
                const arma::vec & cmult,
                const double & tau,
                const double & tau_ext,
                const int & nvar,
                const int & nvar_total,
                const NumericVector & upper_cl,
                const NumericVector & lower_cl,
                const int & ne,
                const int & nx,
                NumericVector & lam_cur,
                NumericVector & lam_prev,
                LogicalVector & strong,
                std::vector<int> & active,
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

    // initialize vars
    int nin = 0;
    bool jerr = 0;
    arma::vec lasso_part(nvar_total);
    arma::vec quantile_part(nvar_total);
    arma::vec ridge_part(nvar_total);

    // compute penalty components and strong rules
    double lam_diff = 2.0 * lam_cur[0] - lam_prev[0];
    double qnt = 2 * tau - 1;
    for (int k = 0; k < nvar; ++k) {
        lasso_part[k] = ptype[k] * lam_cur[0];
        quantile_part[k] = qnt;
        ridge_part[k] = xv[k] + (cmult[k] - ptype[k]) * lam_cur[0];
        if (!strong[k]) {
            strong[k] = std::abs(g[k]) > ptype[k] * lam_diff * std::abs(qnt + (g[k] > 0.0 ? 1.0 : -1.0));
        }
    }

    lam_diff = 2.0 * lam_cur[1] - lam_prev[1];
    qnt = 2 * tau_ext - 1;
    for (int k = nvar; k < nvar_total; ++k) {
        lasso_part[k] = ptype[k] * lam_cur[1];
        quantile_part[k] = qnt;
        ridge_part[k] = xv[k] + (cmult[k] - ptype[k]) * lam_cur[1];
        if (!strong[k]) {
            strong[k] = std::abs(g[k]) > ptype[k] * lam_diff * std::abs(qnt + (g[k] > 0.0 ? 1.0 : -1.0));
        }
    }

    bool kkt_satisfied = false;
    while(!kkt_satisfied && !jerr) {
        bool converge_final = false;
        while(!converge_final && !jerr) {
            double dlx = 0.0;
            for (int k = 0; k < nvar_total; ++k) {
                if (strong[k]) {
                    double gk = arma::dot(x.unsafe_col(k), resid % w);
                    double bk = b[k];
                    double u = gk + bk * xv[k];
                    double v = std::abs(u) - lasso_part[k] * std::abs(quantile_part[k] + (u > 0.0 ? 1.0 : -1.0));
                    if (v > 0.0) {
                        b[k] = std::max(lower_cl[k], std::min(upper_cl[k], copysign(v, u) / ridge_part[k]));
                    } else {
                        b[k] = 0.0;
                    }
                    if (b[k] != bk) {
                        if (mm[k] == 0) {
                            active.push_back(k);
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
                }
            }
            nlp += 1;
            if (nlp > maxit) {
                jerr = 1;
                errcode = -10000;
            }
            if (dlx < thr) {
                converge_final = true;
                // check kkt conditions
                kkt_satisfied = true;
                for (int k = 0; k < nvar_total; k++) {
                    if (!strong[k]) {
                        g[k] = arma::dot(x.unsafe_col(k), resid % w);
                        if (std::abs(g[k]) > lasso_part[k] * std::abs(quantile_part[k] + (g[k] > 0.0 ? 1.0 : -1.0))) {
                            strong[k] = true;
                            kkt_satisfied = false;
                        }
                    }
                }
            } else {
                bool converge = false;
                while(!converge && !jerr) {
                    double dlx = 0.0;
                    for (std::vector<int>::const_iterator k = active.begin(); k != active.end(); ++k) {
                        double gk = arma::dot(x.unsafe_col(*k), resid % w);
                        double bk = b[*k];
                        double u = gk + bk * xv[*k];
                        double v = std::abs(u) - lasso_part[*k] * std::abs(quantile_part[*k] + (u > 0.0 ? 1.0 : -1.0));
                        if (v > 0.0) {
                            b[*k] = std::max(lower_cl[*k], std::min(upper_cl[*k], copysign(v, u) / ridge_part[*k]));
                        } else {
                            b[*k] = 0.0;
                        }
                        if (b[*k] != bk) {
                            double del = b[*k] - bk;
                            resid -= del * x.unsafe_col(*k);
                            dev_cur += del * (2.0 * gk - del * xv[*k]);
                            dlx = std::max(xv[*k] * del * del, dlx);
                        }
                    }
                    nlp += 1;
                    if (nlp > maxit) {
                        jerr = 1;
                        errcode = -10000;
                    }
                    if (dlx < thr) {
                        converge = true;
                    }
                }
            }
        }
    }
    coef.col(idx_lam) = b;
    dev[idx_lam] = dev_cur;
}
