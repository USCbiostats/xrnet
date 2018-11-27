#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#include "hierr_utils.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * Coordinate descent function to fit model
 * for design matrix x on y given
 * variable-specific penalty values and penalty types
 */

// [[Rcpp::export]]
void coord_desc(const arma::mat & x,
                arma::vec & resid,
                const arma::vec & w,
                const arma::vec & ptype,
                const arma::vec & cmult,
                const NumericVector & qnt,
                const NumericVector & lam_cur,
                const IntegerVector & stop,
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
                IntegerVector & nin) {

    bool jerr = 0;
    const int nblocks = lam_cur.length();
    bool kkt_satisfied = false;
    while(!kkt_satisfied && !jerr) {
        bool converge_final = false;
        while(!converge_final && !jerr) {
            double dlx = 0.0;
            int begin = 0;
            for (int blk = 0; blk < nblocks; ++blk) {
                double lambda = lam_cur[blk];
                double q = qnt[blk];
                double end = stop[blk];
                for (int k = begin; k < end; ++k) {
                    if (strong[k]) {
                        double gk = arma::dot(x.unsafe_col(k), resid % w);
                        double bk = b[k];
                        double u = gk + bk * xv[k];
                        double v = std::abs(u) - ptype[k] * cmult[k] * lambda * std::abs(q + sgn(u));
                        if (v > 0.0) {
                            b[k] = std::max(lcl[k],
                                        std::min(ucl[k],
                                        copysign(v, u) / (xv[k] + cmult[k] * (1 - ptype[k]) * lambda)));
                        } else {
                            b[k] = 0.0;
                        }
                        if (b[k] != bk) {
                            if (!active[k]) {
                                ++nin[blk];
                                active[k] = true;
                            }
                            if (sum(nin) <= nx) {
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
                begin = end;
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
                int begin = 0;
                for (int blk = 0; blk < nblocks; ++blk) {
                    double lambda = lam_cur[blk];
                    double q = qnt[blk];
                    double end = stop[blk];
                    for (int k = begin; k < end; ++k) {
                        if (!strong[k]) {
                            g[k] = arma::dot(x.unsafe_col(k), resid % w);
                            if (std::abs(g[k]) > ptype[k] * cmult[k] * lambda * std::abs(q + sgn(g[k]))) {
                                strong[k] = true;
                                kkt_satisfied = false;
                            }
                        }
                    }
                    begin = end;
                }
            } else {
                bool converge = false;
                while(!converge && !jerr) {
                    double dlx = 0.0;
                    int begin = 0;
                    for (int blk = 0; blk < nblocks; ++blk) {
                        double lambda = lam_cur[blk];
                        double q = qnt[blk];
                        double end = stop[blk];
                        for (int k = begin; k < end; ++k) {
                            if (active[k]) {
                                double gk = arma::dot(x.unsafe_col(k), resid  % w);
                                double bk = b[k];
                                double u = gk + bk * xv[k];
                                double v = std::abs(u) - ptype[k] * cmult[k] * lambda * std::abs(q + sgn(u));
                                if (v > 0.0) {
                                    b[k] = std::max(lcl[k],
                                                std::min(ucl[k],
                                                copysign(v, u) / (xv[k] + cmult[k] * (1 - ptype[k]) * lambda)));
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
                        begin = end;
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
