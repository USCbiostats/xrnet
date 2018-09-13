#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#include "hierr_utils.h"

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
                const arma::vec & lasso_part,
                const arma::vec & ridge_part,
                const double & qx,
                const double & qext,
                const int & nv_x,
                const int & nvar_total,
                const NumericVector & ucl,
                const NumericVector & lcl,
                const int & ne,
                const int & nx,
                LogicalVector & strong,
                LogicalVector & active,
                const double & thr,
                const int & maxit,
                const arma::vec & xv,
                arma::vec & b,
                arma::vec & g,
                double & dev_cur,
                double & errcode,
                int & nlp,
                int & nin_x,
                int & nin_ext) {

    bool jerr = 0;
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
                    double v = std::abs(u) - lasso_part[k] * std::abs(qx + sgn(u));
                    if (v > 0.0) {
                        b[k] = std::max(lcl[k], std::min(ucl[k], copysign(v, u) * ridge_part[k]));
                    } else {
                        b[k] = 0.0;
                    }
                    if (b[k] != bk) {
                        if (!active[k]) {
                            //active_x[nin_x] = k;
                            ++nin_x;
                            active[k] = true;
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
                    double v = std::abs(u) - lasso_part[k] * std::abs(qext + sgn(u));
                    if (v > 0.0) {
                        b[k] = std::max(lcl[k], std::min(ucl[k], copysign(v, u) * ridge_part[k]));
                    } else {
                        b[k] = 0.0;
                    }
                    if (b[k] != bk) {
                        if (!active[k]) {
                            //active_ext[nin_ext] = k;
                            ++nin_ext;
                            active[k] = true;
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
                for (int k = 0; k < nv_x; ++k) {
                    if (!strong[k]) {
                        g[k] = arma::dot(x.unsafe_col(k), resid % w);
                        if (std::abs(g[k]) > lasso_part[k] * std::abs(qx + sgn(g[k]))) {
                            strong[k] = true;
                            kkt_satisfied = false;
                        }
                    }
                }
                for (int k = nv_x; k < nvar_total; ++k) {
                    if (!strong[k]) {
                        g[k] = arma::dot(x.unsafe_col(k), resid % w);
                        if (std::abs(g[k]) > lasso_part[k] * std::abs(qext + sgn(g[k]))) {
                            strong[k] = true;
                            kkt_satisfied = false;
                        }
                    }
                }
            } else {
                bool converge = false;
                while(!converge && !jerr) {
                    double dlx = 0.0;
                    /*
                    IntegerVector::const_iterator it;
                    for (it = active_x.begin(); it != active_x.begin() + nin_x; ++it) {
                        int cx = *it;
                        double gk = arma::dot(x.unsafe_col(cx), resid % w);
                        double bk = b[cx];
                        double u = gk + bk * xv[cx];
                        double v = std::abs(u) - lasso_part[cx] * std::abs(qx + sgn(u));
                        if (v > 0.0) {
                            b[cx] = std::max(lcl[cx], std::min(ucl[cx], copysign(v, u) / ridge_part[cx]));
                        } else {
                            b[cx] = 0.0;
                        }
                        if (b[cx] != bk) {
                            double del = b[cx] - bk;
                            resid -= del * x.unsafe_col(cx);
                            dev_cur += del * (2.0 * gk - del * xv[cx]);
                            dlx = std::max(xv[cx] * del * del, dlx);
                        }
                    }
                    for (it = active_ext.begin(); it != active_ext.begin() + nin_ext; ++it) {
                        int cx = *it;
                        double gk = arma::dot(x.unsafe_col(cx), resid % w);
                        double bk = b[cx];
                        double u = gk + bk * xv[cx];
                        double v = std::abs(u) - lasso_part[cx] * std::abs(qext + sgn(u));
                        if (v > 0.0) {
                            b[cx] = std::max(lcl[cx], std::min(ucl[cx], copysign(v, u) / ridge_part[cx]));
                        } else {
                            b[cx] = 0.0;
                        }
                        if (b[cx] != bk) {
                            double del = b[cx] - bk;
                            resid -= del * x.unsafe_col(cx);
                            dev_cur += del * (2.0 * gk - del * xv[cx]);
                            dlx = std::max(xv[cx] * del * del, dlx);
                        }
                    }
                    */
                    for (int k = 0; k < nv_x; ++k) {
                        if (active[k]) {
                            double gk = arma::dot(x.unsafe_col(k), resid % w);
                            double bk = b[k];
                            double u = gk + bk * xv[k];
                            double v = std::abs(u) - lasso_part[k] * std::abs(qx + sgn(u));
                            if (v > 0.0) {
                                b[k] = std::max(lcl[k], std::min(ucl[k], copysign(v, u) * ridge_part[k]));
                            } else {
                                b[k] = 0.0;
                            }
                            if (b[k] != bk) {
                                double del = b[k] - bk;
                                resid -= del * x.unsafe_col(k);
                                dev_cur += del * (2.0 * gk - del * xv[k]);
                                dlx = std::max(xv[k] * del * del, dlx);
                            }
                        }
                    }
                    for (int k = nv_x; k < nvar_total; ++k) {
                        if (active[k]) {
                            double gk = arma::dot(x.unsafe_col(k), resid % w);
                            double bk = b[k];
                            double u = gk + bk * xv[k];
                            double v = std::abs(u) - lasso_part[k] * std::abs(qext + sgn(u));
                            if (v > 0.0) {
                                b[k] = std::max(lcl[k], std::min(ucl[k], copysign(v, u) * ridge_part[k]));
                            } else {
                                b[k] = 0.0;
                            }
                            if (b[k] != bk) {
                                double del = b[k] - bk;
                                resid -= del * x.unsafe_col(k);
                                dev_cur += del * (2.0 * gk - del * xv[k]);
                                dlx = std::max(xv[k] * del * del, dlx);
                            }
                        }
                    }
                    ++nlp;
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
}
