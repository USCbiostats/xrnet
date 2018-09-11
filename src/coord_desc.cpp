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
                const arma::vec & ptype_ind,
                const arma::vec & cmult,
                const double & qx,
                const double & qext,
                const int & nv_x,
                const int & nvar_total,
                const NumericVector & upper_cl,
                const NumericVector & lower_cl,
                const int & ne,
                const int & nx,
                NumericVector & lam_cur,
                NumericVector & lam_prev,
                LogicalVector & strong,
                IntegerVector & active_x,
                IntegerVector & active_ext,
                const double & thr,
                const int & maxit,
                const arma::vec & xv,
                arma::mat & coef,
                arma::vec & b,
                arma::vec & g,
                NumericVector & dev,
                double & dev_cur,
                LogicalVector & ever_active,
                double & errcode,
                int & nlp,
                int & idx_lam,
                int & nin_x,
                int & nin_ext) {

    // initialize vars
    bool jerr = 0;
    arma::vec lasso_part(nvar_total);
    arma::vec ridge_part(nvar_total);

    // compute penalty components and strong rules
    double lam_diff = 2.0 * lam_cur[0] - lam_prev[0];
    for (int k = 0; k < nv_x; ++k) {
        lasso_part[k] = ptype_ind[k] * lam_cur[0];
        ridge_part[k] = xv[k] + (cmult[k] - ptype_ind[k]) * lam_cur[0];
        if (!strong[k]) {
            strong[k] = std::abs(g[k]) > ptype_ind[k] * lam_diff * std::abs(qx + (g[k] > 0.0 ? 1.0 : -1.0));
        }
    }

    lam_diff = 2.0 * lam_cur[1] - lam_prev[1];
    for (int k = nv_x; k < nvar_total; ++k) {
        lasso_part[k] = ptype_ind[k] * lam_cur[1];
        ridge_part[k] = xv[k] + (cmult[k] - ptype_ind[k]) * lam_cur[1];
        if (!strong[k]) {
            strong[k] = std::abs(g[k]) > ptype_ind[k] * lam_diff * std::abs(qext + (g[k] > 0.0 ? 1.0 : -1.0));
        }
    }

    bool kkt_satisfied = false;
    while(!kkt_satisfied && !jerr) {
        bool converge_final = false;
        while(!converge_final && !jerr) {
            double dlx = 0.0;
            for (int k = 0; k < nv_x; ++k) {
                if (strong[k]) {
                    double gk = arma::dot(x.unsafe_col(k), resid % w);
                    double bk = b[k];
                    double u = gk + bk * xv[k];
                    double v = std::abs(u) - lasso_part[k] * std::abs(qx + (u > 0.0 ? 1.0 : -1.0));
                    if (v > 0.0) {
                        b[k] = std::max(lower_cl[k], std::min(upper_cl[k], copysign(v, u) / ridge_part[k]));
                    } else {
                        b[k] = 0.0;
                    }
                    if (b[k] != bk) {
                        if (!ever_active[k]) {
                            active_x[nin_x] = k;
                            ++nin_x;
                            ever_active[k] = true;
                        }
                        if (nin_x + nin_ext <= nx) {
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
            for (int k = nv_x; k < nvar_total; ++k) {
                if (strong[k]) {
                    double gk = arma::dot(x.unsafe_col(k), resid % w);
                    double bk = b[k];
                    double u = gk + bk * xv[k];
                    double v = std::abs(u) - lasso_part[k] * std::abs(qext + (u > 0.0 ? 1.0 : -1.0));
                    if (v > 0.0) {
                        b[k] = std::max(lower_cl[k], std::min(upper_cl[k], copysign(v, u) / ridge_part[k]));
                    } else {
                        b[k] = 0.0;
                    }
                    if (b[k] != bk) {
                        if (!ever_active[k]) {
                            active_ext[nin_ext] = k;
                            ++nin_ext;
                            ever_active[k] = true;
                        }
                        if (nin_x + nin_ext <= nx) {
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
            ++nlp;
            if (nlp > maxit) {
                jerr = 1;
                errcode = -10000;
            }
            if (dlx < thr) {
                converge_final = true;
                // check kkt conditions
                kkt_satisfied = true;
                for (int k = 0; k < nv_x; k++) {
                    if (!strong[k]) {
                        g[k] = arma::dot(x.unsafe_col(k), resid % w);
                        if (std::abs(g[k]) > lasso_part[k] * std::abs(qx + (g[k] > 0.0 ? 1.0 : -1.0))) {
                            strong[k] = true;
                            kkt_satisfied = false;
                        }
                    }
                }
                for (int k = nv_x; k < nvar_total; k++) {
                    if (!strong[k]) {
                        g[k] = arma::dot(x.unsafe_col(k), resid % w);
                        if (std::abs(g[k]) > lasso_part[k] * std::abs(qext + (g[k] > 0.0 ? 1.0 : -1.0))) {
                            strong[k] = true;
                            kkt_satisfied = false;
                        }
                    }
                }
            } else {
                bool converge = false;
                while(!converge && !jerr) {
                    double dlx = 0.0;
                    for (IntegerVector::const_iterator k = active_x.begin(); k != active_x.begin() + nin_x; ++k) {
                        double gk = arma::dot(x.unsafe_col(*k), resid % w);
                        double bk = b[*k];
                        double u = gk + bk * xv[*k];
                        double v = std::abs(u) - lasso_part[*k] * std::abs(qx + (u > 0.0 ? 1.0 : -1.0));
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
                    for (IntegerVector::const_iterator k = active_ext.begin(); k != active_ext.begin() + nin_ext; ++k) {
                        double gk = arma::dot(x.unsafe_col(*k), resid % w);
                        double bk = b[*k];
                        double u = gk + bk * xv[*k];
                        double v = std::abs(u) - lasso_part[*k] * std::abs(qext + (u > 0.0 ? 1.0 : -1.0));
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
