#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
#include "create_data.h"
#include "coord_desc.h"

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
                max_pen = std::max(max_pen, g[j] / cmult[j]);
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

    for (int j = 0; j < nlam_total; ++j) {
        for (int i = 0; i < nvar; ++i) {
            double z_alpha = a0[j];
            for (int k = ext_start; k < nvar_total; ++k) {
                z_alpha += coef.at(k, j) * (ext_.at(i, k - ext_start) - xm[k]) / xs[k];
            }
            coef.at(i, j) = ys * (z_alpha + coef.at(i, j)) / xs[i];
        }
    }
}

/*
 * Main function to fit hierarchical regularized
 * regression for gaussian outcome
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
                  const NumericVector & ptype,
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
                  const bool & isd,
                  const bool & isd_ext,
                  const bool & intr,
                  const bool & intr_ext) {

    // Create single level regression matrix
    int nvar_total = nvar + nvar_ext + nvar_unpen + intr_ext;
    arma::vec xm(nvar_total, arma::fill::zeros);
    arma::vec xv(nvar_total, arma::fill::ones);
    arma::vec xs(nvar_total, arma::fill::ones);
    const arma::vec wgt = w / sum(w);
    int ext_start;
    const arma::mat xnew = create_data(nobs, nvar, nvar_ext, nvar_unpen, nvar_total, x_, ext_, fixed_,
                                       wgt, isd, isd_ext, intr, intr_ext, xm, xv, xs, ext_start);

    // determine non-constant variables -- still to be done

    // standardize y, confidence limits, user penalties
    double ym = 0.0;
    double ys = 1.0;
    arma::vec outer_resid = y_;
    standardize_vec(outer_resid, wgt, ym, ys, intr);
    lower_cl = lower_cl / ys;
    upper_cl = upper_cl / ys;
    NumericVector ulam = Rcpp::clone(ulam_);
    NumericVector ulam_ext = Rcpp::clone(ulam_ext_);
    //ulam = ulam_ / ys;
    //ulam_ext = ulam_ext_ / ys;

    if (isd) {
        for (int i = 0; i < nvar; ++i) {
            if (lower_cl[i] != R_NegInf) {
                lower_cl[i] *= xs[i];
            }
            if (upper_cl[i] != R_PosInf) {
                upper_cl[i] *= upper_cl[i] * xs[i];
            }
        }
    }
    if (isd_ext) {
        for (int i = ext_start; i < nvar_total; ++i) {
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

    // compute gradient vector
    arma::vec g(nvar_total, arma::fill::zeros);
    g = arma::abs(xnew.t() * (wgt % outer_resid));

    // compute penalty paths
    int start = 0;
    NumericVector lam_path = compute_penalty(ulam, nlam, ptype[0], pratio, g, cmult, start, nvar);
    NumericVector lam_path_ext = 0.0;
    if (nvar_ext > 0) {
        lam_path_ext = compute_penalty(ulam_ext, nlam_ext, ptype[ext_start],
                                       pratio_ext, g, cmult, ext_start, nvar_total);
    }

    // loop through all penalty combinations
    NumericVector lam_cur(2, 0.0);
    double dev_outer = 0.0;
    double dev_inner = 0.0;
    double dev_old = 0.0;
    arma::vec inner_resid(nobs);
    arma::vec coef_outer(nvar_total, arma::fill::zeros);
    arma::vec coef_inner(nvar_total, arma::fill::zeros);
    IntegerVector mm(nvar_total, 0);

    int idx_lam = 0;
    int nlp = 0;
    double errcode = 0.0;

    for (int m = 0; m < nlam; ++m) {
        lam_cur[0] = lam_path[m];

        for (int m2 = 0; m2 < nlam_ext; ++m2) {
            lam_cur[1] = lam_path_ext[m2];

            if (m2 == 0) {
                coord_desc(xnew, outer_resid, wgt, ptype, tau, tau_ext, nobs, nvar, nvar_total,
                           cmult, upper_cl, lower_cl, ne, nx,
                           lam_cur, thr, maxit, xv, coef, coef_outer, g,
                           dev, dev_outer, mm, errcode, nlp, idx_lam);

                inner_resid = outer_resid;
                coef_inner = coef_outer;
                dev_inner = dev_outer;
                dev_old = dev_outer;

            }
            else {
                coord_desc(xnew, inner_resid, wgt, ptype, tau, tau_ext, nobs, nvar, nvar_total,
                           cmult, upper_cl, lower_cl, ne, nx,
                           lam_cur, thr, maxit, xv, coef, coef_inner, g,
                           dev, dev_inner, mm, errcode, nlp, idx_lam);

                //stop if max deviance or no appreciable change in deviance
                if (pratio_ext > 0.0) {
                    if ((dev_inner - dev_old) < (1e-05 * dev_inner) || dev_inner > 0.999 || errcode > 0.0) {
                        idx_lam += nlam_ext - m2;
                        break;
                    }
                    else {
                        dev_old = dev_inner;
                    }
                }
            }
            ++idx_lam;
        }
        if (errcode > 0.0) {
            break;
        }
    }

    // compute and unstandardize predictor variables
    arma::vec b0(nlam_total, arma::fill::zeros);
    arma::vec a0(nlam_total, arma::fill::zeros);
    if (intr_ext) {
        a0 = arma::conv_to<arma::colvec>::from(coef.row(nvar + nvar_unpen));
    }
    compute_coef(coef, ext_, nvar, nvar_ext, nvar_total, nlam_total, xm, xs, ys, a0, intr_ext, ext_start);

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

    // fix first penalties
    if (ulam_[0] == 0.0) {
        lam_path[0] = exp(2 * log(lam_path[1]) - log(lam_path[2]));
    }
    if (ulam_ext_[0] == 0.0) {
        lam_path_ext[0] = exp(2 * log(lam_path_ext[1]) - log(lam_path_ext[2]));
    }

    // return model fit for all penalty combinations
    return Rcpp::List::create(Named("beta0") = b0,
                              Named("betas") = coef.head_rows(nvar),
                              Named("gammas") = gammas,
                              Named("alpha0") = a0,
                              Named("alphas") = alphas,
                              Named("penalty") = lam_path,
                              Named("penalty_ext") = lam_path_ext,
                              Named("penalty_type") = ptype[0],
                              Named("quantile") = tau,
                              Named("penalty_type_ext") = ptype[nvar_total - 1],
                              Named("quantile_ext") = tau_ext,
                              Named("penalty_ratio") = pratio,
                              Named("penalty_ratio_ext") = pratio_ext,
                              Named("deviance") = dev,
                              Named("nlp") = nlp);
}
