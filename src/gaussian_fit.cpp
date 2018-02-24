#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * Generates matrix [x | x*external]
 */

// [[Rcpp::export]]
arma::mat create_data(const int & nobs,
                      const int & nvar,
                      const int & nvar_ext,
                      const int & nvar_unpen,
                      const int & nvar_total,
                      const arma::mat & x,
                      const arma::mat & ext,
                      const arma::mat & unpen,
                      const arma::vec & w,
                      const bool & isd,
                      const bool & isd_ext,
                      const bool & intr,
                      const bool & intr_ext,
                      arma::vec & xm,
                      arma::vec & xv,
                      arma::vec & xs,
                      int & ext_start) {

    // initialize new design matrix
    arma::mat xnew(nobs, nvar_total);

    // standardize x (predictor variables)
    if (intr) {
        if (isd) {
            for (int j = 0; j < nvar; ++j) {
                xm[j] = arma::dot(w, x.unsafe_col(j));
                xs[j] = sqrt(arma::dot(x.unsafe_col(j) - xm[j], (x.unsafe_col(j) - xm[j]) % w));
                xnew.col(j) = (x.unsafe_col(j) - xm[j]) / xs[j];
            }
        }
        else {
            for (int j = 0; j < nvar; ++j) {
                xm[j] = arma::dot(w, x.unsafe_col(j));
                xnew.col(j) = x.unsafe_col(j) - xm[j];
                xv[j] = arma::dot(xnew.unsafe_col(j), xnew.unsafe_col(j) % w);
            }
        }
    }
    else {
        if (isd) {
            for (int j = 0; j < nvar; ++j) {
                double xm_j = arma::dot(w, x.unsafe_col(j));
                double vc = arma::dot(x.unsafe_col(j) - xm_j, (x.unsafe_col(j) - xm_j) % w);
                xs[j] = sqrt(vc);
                xnew.col(j) = x.unsafe_col(j) / xs[j];
                xv[j] = 1.0 + xm_j * xm_j / vc;
            }
        }
        else {
            for (int j = 0; j < nvar; ++j) {
                xnew.col(j) = x.unsafe_col(j);
                xv[j] = arma::dot(x.unsafe_col(j), x.unsafe_col(j) % w);
            }
        }
    }

    // Create reference to standardized x variables in xnew
    arma::mat xsub(xnew.memptr(), nobs, nvar, false, true);

    // standardize unpenalized variables (unpen)
    int xnew_col = nvar;
    if (nvar_unpen > 0) {
        if (intr) {
            if (isd) {
                for (int j = 0; j < nvar_unpen; ++j) {
                    xm[xnew_col] = arma::dot(w, unpen.unsafe_col(j));
                    xs[xnew_col] = sqrt(arma::dot(unpen.unsafe_col(j) - xm[xnew_col], (unpen.unsafe_col(j) - xm[xnew_col]) % w));
                    xnew.col(xnew_col) = (unpen.unsafe_col(j) - xm[xnew_col]) / xs[xnew_col];
                    ++xnew_col;
                }
            }
            else {
                for (int j = 0; j < nvar_unpen; ++j) {
                    xm[xnew_col] = arma::dot(w, unpen.unsafe_col(j));
                    xnew.col(xnew_col) = unpen.unsafe_col(j) - xm[xnew_col];
                    xv[xnew_col] = arma::dot(xnew.unsafe_col(xnew_col), xnew.unsafe_col(xnew_col) % w);
                    ++xnew_col;
                }
            }

        }
        else {
            if (isd) {
                for (int j = 0; j < nvar_unpen; ++j) {
                    double xm_j = arma::dot(w, unpen.unsafe_col(j));
                    double vc = arma::dot(unpen.unsafe_col(j) - xm_j, (unpen.unsafe_col(j) - xm_j) % w);
                    xs[xnew_col] = sqrt(vc);
                    xnew.col(xnew_col) = unpen.unsafe_col(j) / xs[xnew_col];
                    xv[xnew_col] = 1.0 + xm_j * xm_j / vc;
                    ++xnew_col;
                }
            }
            else {
                for (int j = 0; j < nvar_unpen; ++j) {
                    xnew.col(xnew_col) = unpen.unsafe_col(j);
                    xv[xnew_col] = arma::dot(unpen.unsafe_col(j), unpen.unsafe_col(j) % w);
                    ++xnew_col;
                }
            }

        }
    }

    // add 2nd level intercept column
    if (intr_ext) {
        xnew.col(xnew_col) = xsub * arma::ones<arma::mat>(nvar, 1);
        xv[xnew_col] = arma::dot(xnew.unsafe_col(xnew_col), xnew.unsafe_col(xnew_col)) / nobs;
        xs[xnew_col] = 1.0;
        ++xnew_col;
    }

    // Standardize external variables (ext)
    ext_start = xnew_col;
    if (nvar_ext > 0) {
        if (intr_ext) {
            if (isd_ext) {
                for (int j = 0; j < nvar_ext; ++j) {
                    xm[xnew_col] = arma::mean(ext.unsafe_col(j));
                    xs[xnew_col] = sqrt(arma::dot(ext.unsafe_col(j) - xm[xnew_col], ext.unsafe_col(j) - xm[xnew_col]) / nvar);
                    xnew.col(xnew_col) = xsub * ((ext.unsafe_col(j) - xm[xnew_col]) / xs[xnew_col]);
                    xv[xnew_col] = arma::dot(xnew.unsafe_col(xnew_col), xnew.unsafe_col(xnew_col)) / nobs;
                    ++xnew_col;
                }
            }
            else {
                for (int j = 0; j < nvar_ext; ++j) {
                    xm[xnew_col] = arma::mean(ext.unsafe_col(j));
                    xnew.col(xnew_col) = xsub * (ext.unsafe_col(j) - xm[xnew_col]);
                    xv[xnew_col] = arma::dot(xnew.unsafe_col(xnew_col), xnew.unsafe_col(xnew_col)) / nobs;
                    ++xnew_col;
                }
            }
        }
        else {
            if (isd_ext) {
                for (int j = 0; j < nvar_ext; ++j) {
                    double xm_j = arma::mean(ext.unsafe_col(j));
                    xs[xnew_col] = sqrt(arma::dot(ext.unsafe_col(j) - xm_j, ext.unsafe_col(j) - xm_j) / nvar);
                    xnew.col(xnew_col) = xsub * (ext.unsafe_col(j) / xs[xnew_col]);
                    xv[xnew_col] = arma::dot(xnew.unsafe_col(xnew_col), xnew.unsafe_col(xnew_col)) / nobs;
                    ++xnew_col;
                }
            }
            else {
                for (int j = 0; j < nvar_ext; ++j) {
                    xnew.col(xnew_col) = xsub * ext.unsafe_col(j);
                    xv[xnew_col] = arma::dot(xnew.unsafe_col(xnew_col), xnew.unsafe_col(xnew_col)) / nobs;
                    ++xnew_col;
                }
            }
        }
    }
    return(xnew);
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

    int nvar_total = nvar + nvar_ext + nvar_unpen;
    if (intr_ext) {++nvar_total;}

    // Create single level regression matrix
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

    // initialize objects to hold fitting results
    int nlam_total = nlam * nlam_ext;
    arma::mat coef(nvar_total, nlam_total, arma::fill::zeros);
    NumericVector dev(nlam_total);

    // ---------run coordinate descent for all penalties ----------

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
