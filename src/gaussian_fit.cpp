#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#include "hierr_utils.h"
#include "create_data.h"
#include "coord_desc.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * Main function to fit hierarchical regularized
 * regression for gaussian outcome (dense external data)
 */

// [[Rcpp::export]]
List gaussian_fit(const arma::mat & x_,
                  const arma::vec & y_,
                  const arma::mat & ext_,
                  const arma::mat & fixed_,
                  const int & nobs,
                  const int & nvar,
                  const int & nvar_ext,
                  const int & nvar_unpen,
                  const arma::vec & w,
                  const arma::vec & ptype,
                  const double & tau,
                  const double & tau_ext,
                  const arma::vec & cmult,
                  NumericVector & lower_cl,
                  NumericVector & upper_cl,
                  const int & ne,
                  const int & nx,
                  const int & nlam,
                  const int & nlam_ext,
                  const double & pratio,
                  const double & pratio_ext,
                  NumericVector ulam_,
                  NumericVector ulam_ext_,
                  const double & thr,
                  const int & maxit,
                  const bool & earlyStop,
                  const bool & isd,
                  const bool & isd_ext,
                  const bool & intr,
                  const bool & intr_ext) {

    // Create single level regression matrix
    const int ext_start = nvar + nvar_unpen + intr_ext;
    const int nvar_total = ext_start + nvar_ext;
    const int nv_x = nvar + nvar_unpen;
    const int nv_ext = intr_ext + nvar_ext;
    arma::vec xm(nvar_total, arma::fill::zeros);
    arma::vec xv(nvar_total, arma::fill::ones);
    arma::vec xs(nvar_total, arma::fill::ones);
    const arma::vec wgt = w / sum(w);
    const arma::mat xnew = create_data(nobs, nvar, nvar_ext, nvar_unpen,
                                       nvar_total, x_, ext_, fixed_,
                                       wgt, isd, isd_ext, intr, intr_ext,
                                       xm, xv, xs);

    // determine non-constant variables -- still to be done

    // standardize y, confidence limits, user penalties
    double ym = 0.0;
    double ys = 1.0;
    arma::vec outer_resid = y_;
    standardize_vec(outer_resid, wgt, ym, ys, intr);
    lower_cl = lower_cl / ys;
    upper_cl = upper_cl / ys;
    NumericVector ulam = Rcpp::clone(ulam_) / ys;
    NumericVector ulam_ext = Rcpp::clone(ulam_ext_) / ys;

    if (isd) {
        for (int i = 0; i < nv_x; ++i) {
            if (lower_cl[i] != R_NegInf) {
                lower_cl[i] *= xs[i];
            }
            if (upper_cl[i] != R_PosInf) {
                upper_cl[i] *= upper_cl[i] * xs[i];
            }
        }
    }
    if (isd_ext && nv_ext > 0) {
        for (int i = nv_x; i < nvar_total; ++i) {
            if (lower_cl[i] != R_NegInf) {
                lower_cl[i] *= xs[i];
            }
            if (upper_cl[i] != R_PosInf) {
                upper_cl[i] *= xs[i];
            }
        }
    }

    // ---------run coordinate descent for all penalties ----------

    // initialize objects to hold fitting results
    int nlam_total = nlam * nlam_ext;
    arma::mat coef(nvar_total, nlam_total, arma::fill::zeros);
    NumericVector dev(nlam_total);
    NumericVector num_passes(nlam_total);

    // compute gradient vector
    arma::vec gouter = xnew.t() * (wgt % outer_resid);

    // compute penalty paths
    NumericVector lam_path(nlam);
    compute_penalty(lam_path, ulam, ptype[0],
                    pratio, gouter, cmult, 0, nvar);

    NumericVector lam_path_ext(nlam_ext);
    if (nvar_ext > 0) {
        compute_penalty(lam_path_ext, ulam_ext,
                        ptype[ext_start], pratio_ext,
                        gouter, cmult, ext_start, nvar_total);
    } else {
        lam_path_ext[0] = 0.0;
    }

    NumericVector lam_cur(2, 0.0);
    NumericVector lam_prev(2, 0.0);
    double dev_outer = 0.0, dev_inner = 0.0, dev_old = 0.0;
    double errcode = 0.0;
    int nlp_old = 0, idx_lam = 0, nlp = 0;
    arma::vec inner_resid(nobs);
    arma::vec coef_outer(nvar_total, arma::fill::zeros);
    arma::vec coef_inner(nvar_total, arma::fill::zeros);
    arma::vec ginner(nvar_total);

    // vars to track strong set and active set
    LogicalVector ever_active(nvar_total, 0);
    LogicalVector strong(nvar_total, 0);
    IntegerVector nin(2, 0);
    //IntegerVector active_x(nv_x, 0);
    //IntegerVector active_ext(nv_ext, 0);

    // quantile constants
    const NumericVector qnt = NumericVector::create(2 * tau - 1, 2 * tau_ext - 1);
    const IntegerVector blkend = IntegerVector::create(nv_x, nvar_total);

    // loop through all penalty combinations
    for (int m = 0; m < nlam; ++m) {
        lam_cur[0] = lam_path[m];

        for (int m2 = 0; m2 < nlam_ext; ++m2) {
            lam_cur[1] = lam_path_ext[m2];

            if (m2 == 0) {
                // reset strong / active for ext vars
                if (nv_ext > 0) {
                    std::fill(ever_active.begin() + nv_x, ever_active.end(), false);
                    std::fill(strong.begin() + nv_x, strong.end(), false);
                    nin[1] = 0;
                }

                // update strong
                updateStrong(strong, gouter, ptype, cmult, lam_cur,
                             lam_prev, qnt, blkend);

                // fit model
                coord_desc(xnew, outer_resid, wgt, ptype, cmult, qnt,
                           lam_cur, blkend, upper_cl, lower_cl, ne, nx,
                           strong, ever_active, thr, maxit, xv,
                           coef_outer, gouter, dev_outer, errcode, nlp, nin);

                // copy results for inner loop
                inner_resid = outer_resid;
                coef_inner = coef_outer;
                ginner = gouter;
                dev_inner = dev_outer;
                dev_old = dev_outer;
            }
            else {
                // update strong
                updateStrong(strong, ginner, ptype, cmult,
                             lam_cur, lam_prev, qnt, blkend);

                // fit model
                coord_desc(xnew, inner_resid, wgt, ptype, cmult,
                           qnt, lam_cur, blkend, upper_cl, lower_cl,
                           ne, nx, strong, ever_active, thr, maxit,
                           xv, coef_inner, ginner, dev_inner, errcode,
                           nlp, nin);
            }

            // save results
            coef.col(idx_lam) = coef_inner;
            dev[idx_lam] = dev_inner;
            num_passes[idx_lam] = nlp - nlp_old;
            nlp_old = nlp;

            // check error
            if (errcode != 0.0) {
                break;
            }

            // check stop conditions
            if (pratio_ext > 0.0 && earlyStop && m2 != 0) {
                double dev_diff = dev_inner - dev_old;
                if (dev_diff < (1e-05 * dev_inner) || dev_inner > 0.999) {
                    idx_lam += nlam_ext - m2;
                    break;
                }
            } else {
                dev_old = dev_inner;
            }

            if (lam_cur[1] == 9.9e35) {
                lam_prev[1] = 0.0;
            } else {
                lam_prev[1] = lam_cur[1];
            }
            ++idx_lam;
        }
        if (lam_cur[0] == 9.9e35) {
            lam_prev[0] = 0.0;
        } else {
            lam_prev[0] = lam_cur[0];
        }
    }

    // compute and unstandardize predictor variables
    arma::vec a0(nlam_total, arma::fill::zeros);
    if (intr_ext) {
        a0 = arma::conv_to<arma::colvec>::from(coef.row(nv_x));
    }
    compute_coef(coef, ext_, nvar, nvar_ext,
                 nvar_total, nlam_total, xm, xs,
                 ys, a0, intr_ext, ext_start);

    // unstandardize unpenalized variables
    arma::mat gammas;
    if (nvar_unpen > 0) {
        for (int j = 0; j < nlam_total; ++j) {
            for (int i = nvar; i < nv_x; ++i) {
                coef.at(i, j) = ys * coef.at(i, j) / xs[i];
            }
        }
        gammas = coef.rows(nvar, nv_x - 1);
    } else {
        gammas = 0.0;
    }

    arma::vec b0(nlam_total, arma::fill::zeros);
    if (intr) {
        b0 = ym - ((xm.head(nv_x)).t() * coef.head_rows(nv_x)).t();
    }

    // unstandardize external variables
    arma::mat alphas;
    if (nvar_ext > 0) {
        for (int j = 0; j < nlam_total; ++j) {
            for (int i = ext_start; i < nvar_total; ++i) {
                coef.at(i, j) = ys * coef.at(i, j) / xs[i];
            }
        }
        if (intr_ext) {
            a0 = (arma::mean(coef.head_rows(nvar)) - (xm.tail(nvar_ext)).t() * coef.tail_rows(nvar_ext)).t();
        }
        alphas = coef.tail_rows(nvar_ext);
    } else {
        alphas = 0.0;
    }

    // fix first penalties (when path automatically computed)
    if (ulam_[0] == 0.0) {
        lam_path[0] = exp(2 * log(lam_path[1]) - log(lam_path[2]));
    }
    if (ulam_ext_[0] == 0.0 && nvar_ext > 0) {
        lam_path_ext[0] = exp(2 * log(lam_path_ext[1]) - log(lam_path_ext[2]));
    }

    // return model fit for all penalty combinations
    return Rcpp::List::create(Named("beta0") = b0,
                              Named("betas") = coef.head_rows(nvar),
                              Named("gammas") = gammas,
                              Named("alpha0") = a0,
                              Named("alphas") = alphas,
                              Named("penalty") = lam_path * ys,
                              Named("penalty_ext") = lam_path_ext * ys,
                              Named("penalty_type") = ptype[0],
                              Named("quantile") = tau,
                              Named("penalty_type_ext") = ptype[nvar_total - 1],
                              Named("quantile_ext") = tau_ext,
                              Named("penalty_ratio") = pratio,
                              Named("penalty_ratio_ext") = pratio_ext,
                              Named("deviance") = dev,
                              Named("num_passes") = num_passes,
                              Named("nlp") = nlp,
                              Named("status") = errcode);
}
