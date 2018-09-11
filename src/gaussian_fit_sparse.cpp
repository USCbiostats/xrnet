#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#include "hierr_utils.h"
#include "create_data_sparse.h"
#include "coord_desc.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


/*
* Main function to fit hierarchical regularized
* regression for gaussian outcome
* external data is in csc sparse matrix format
*/

    // [[Rcpp::export]]
List gaussian_fit_sparse(const arma::mat & x_,
                         const arma::vec & y_,
                         const arma::sp_mat & ext_,
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
    const arma::mat xnew = create_data_sparse(nobs, nvar, nvar_ext, nvar_unpen,
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
    if (isd_ext && nvar_ext > 0) {
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
    arma::vec g = xnew.t() * (wgt % outer_resid);

    // compute individual ptype
    arma::vec ptype_ind = cmult % ptype;

    // compute penalty paths
    int start = 0;
    NumericVector lam_path = compute_penalty(ulam, nlam, ptype[0], pratio,
                                             g, cmult, start, nvar);
    NumericVector lam_path_ext = 0.0;
    if (nvar_ext > 0) {
        lam_path_ext = compute_penalty(ulam_ext, nlam_ext,
                                       ptype[ext_start], pratio_ext,
                                       g, cmult, ext_start, nvar_total);
    }

    NumericVector lam_cur(2, 0.0);
    NumericVector lam_prev(2, 0.0);
    double dev_outer = 0.0, dev_inner = 0.0, dev_old = 0.0;
    double errcode = 0.0;
    int nlp_old = 0, idx_lam = 0, nlp = 0;
    arma::vec inner_resid(nobs);
    arma::vec coef_outer(nvar_total, arma::fill::zeros);
    arma::vec coef_inner(nvar_total, arma::fill::zeros);

    // vars to track strong set and active set
    int nin_x = 0, nin_ext = 0;
    LogicalVector ever_active(nvar_total, 0);
    LogicalVector strong(nvar_total, 0);
    IntegerVector active_x(nv_x, 0);
    IntegerVector active_ext(nv_ext, 0);

    // quantile constants
    const double qx = 2 * tau - 1;
    const double qext = 2 * tau_ext - 1;

    // loop through all penalty combinations
    for (int m = 0; m < nlam; ++m) {
        lam_cur[0] = lam_path[m];

        for (int m2 = 0; m2 < nlam_ext; ++m2) {
            lam_cur[1] = lam_path_ext[m2];

            if (m2 == 0) {

                // reset strong and active for ext vars
                if (nv_ext > 0) {
                    std::fill(ever_active.begin() + nv_x, ever_active.end(), false);
                    std::fill(strong.begin() + nv_x, strong.end(), false);
                    nin_ext = 0;
                }

                coord_desc(xnew, outer_resid, wgt, ptype_ind,
                           cmult, qx, qext, nvar,
                           nvar_total, upper_cl, lower_cl,
                           ne, nx, lam_cur, lam_prev,
                           strong, active_x, active_ext, thr, maxit,
                           xv, coef, coef_outer, g,
                           dev, dev_outer, ever_active, errcode,
                           nlp, idx_lam, nin_x, nin_ext);

                inner_resid = outer_resid;
                coef_inner = coef_outer;
                dev_inner = dev_outer;
                dev_old = dev_outer;
                num_passes[idx_lam] = nlp - nlp_old;
                nlp_old = nlp;

            }
            else {
                coord_desc(xnew, inner_resid, wgt, ptype_ind,
                           cmult, qx, qext, nvar,
                           nvar_total, upper_cl, lower_cl,
                           ne, nx, lam_cur, lam_prev,
                           strong, active_x, active_ext, thr, maxit,
                           xv, coef, coef_inner, g,
                           dev, dev_inner, ever_active, errcode,
                           nlp, idx_lam, nin_x, nin_ext);

                num_passes[idx_lam] = nlp - nlp_old;
                nlp_old = nlp;
                //stop if max deviance or no appreciable change in deviance
                if (pratio_ext > 0.0 && earlyStop) {
                    if ((dev_inner - dev_old) < (1e-05 * dev_inner) || dev_inner > 0.999 || errcode > 0.0) {
                        idx_lam += nlam_ext - m2;
                        break;
                    }
                    else {
                        dev_old = dev_inner;
                    }
                }
            }
            lam_prev[1] = lam_cur[1];
            if (errcode > 0.0) {
                break;
            }
            ++idx_lam;
        }
        lam_prev[0] = lam_cur[0];
    }

    // compute and unstandardize predictor variables
    arma::vec b0(nlam_total, arma::fill::zeros);
    arma::vec a0(nlam_total, arma::fill::zeros);
    if (intr_ext) {
        a0 = arma::conv_to<arma::colvec>::from(coef.row(nvar + nvar_unpen));
    }
    compute_coef_sparse(coef, ext_, nvar, nvar_ext, nvar_total,
                        nlam_total, xm, xs, ys, a0, intr_ext, ext_start);

    // unstandardize unpenalized variables
    arma::mat gammas;
    if (nvar_unpen > 0) {
        for (int j = 0; j < nlam_total; ++j) {
            for (int i = nvar; i < (nvar + nvar_unpen); ++i) {
                coef.at(i, j) = ys * coef.at(i, j) / xs[i];
            }
        }
        gammas = coef.rows(nvar, nvar + nvar_unpen - 1);
    } else {
        gammas = 0.0;
    }

    if (intr) {
        b0 = ym - ((xm.head(nvar + nvar_unpen)).t() * coef.head_rows(nvar + nvar_unpen)).t();
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
