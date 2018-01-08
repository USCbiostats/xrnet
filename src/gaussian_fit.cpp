#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericMatrix create_data(int & nobs,
                          int & nvar,
                          int & nvar_total,
                          arma::mat & x,
                          arma::mat & ext,
                          arma::vec & w,
                          bool & isd,
                          bool & isd_ext,
                          bool & intr,
                          bool & intr_ext,
                          arma::vec & xm,
                          arma::vec & xv,
                          arma::vec & xs) {

    arma::vec v = sqrt(w);
    arma::mat xnew(nobs, nvar_total);

    xm.subvec(0, nvar - 1) = x.t() * w;
    if (intr == true) {
        if (isd == true) {
            for (int j = 0; j < nvar; j++) {
                xnew.col(j) = x.col(j) - xm[j];
                xs[j] = sqrt(arma::dot(xnew.col(j), xnew.col(j)) / nobs);
                xnew.col(j) = xnew.col(j) / xs[j];
                xv[j] = 1.0;
            }
        } else {
            for (int j = 0; j < nvar; j++) {
                xnew.col(j) = x.col(j) - xm[j];
                xv[j] = arma::dot(xnew.col(j), xnew.col(j)) / nobs;
                xs[j] = 1.0;
            }
        }
    } else {
        for (int j = 0; j < nvar; j++) {
            xv[j] = arma::dot(x.col(j), x.col(j)) / nobs;
            if (isd == true) {

            }
        }
    }

    double xm_col;
    double xv_col;
    arma::vec ext_temp(nvar);
    arma::mat xsub = xnew.submat(0, 0, nobs - 1, nvar - 1);

    if (intr_ext == true) {
        xnew.col(nvar) = xnew.submat(0, 0, nobs - 1, nvar - 1) * arma::ones<arma::mat>(nvar, 1);
        xm_col = arma::dot(xnew.col(nvar), w) / nobs;
        xv_col = arma::dot(xnew.col(nvar), xnew.col(nvar)) / nobs;
        xv[nvar] = xv_col - xm_col * xm_col;
        xs[nvar] = 1.0;
        xm.subvec(nvar + 1, nvar_total - 1) = arma::mean(ext).t();
        for (int j = nvar + 1; j < nvar_total; j++) {
            ext_temp = ext.col(j - nvar - 1) - xm[j];
            if (isd_ext == true) {
                xs[j] = sqrt(arma::dot(ext_temp, ext_temp) / nvar);
                ext_temp = ext_temp / xs[j];
            } else {
                xs[j] = 1.0;

            }
            xnew.col(j) = xsub * ext_temp;
            xm_col = arma::dot(xnew.col(j), w) / nobs;
            xv_col = arma::dot(xnew.col(j), xnew.col(j)) / nobs;
            xv[j] = xv_col - xm_col * xm_col;
        }
    }

    return(wrap(xnew));
}

// [[Rcpp::export]]
void standardize_mat(int & nobs,
                     int & nvar,
                     NumericMatrix & x,
                     NumericVector & w,
                     bool & isd,
                     bool & intr,
                     NumericVector & xm,
                     NumericVector & xv,
                     NumericVector & xs) {

    NumericVector v = sqrt(w);

    if (intr == true) {
        for (int j = 0; j < nvar; j++) {
            xm[j] = std::inner_product(x(_, j).begin(), x(_, j).end(), w.begin(), 0.0);
            x(_, j) = x(_, j) - xm[j];
            xv[j] = std::inner_product(x(_, j).begin(), x(_, j).end(), x(_, j).begin(), 0.0) / nobs;
            if (isd) {
                xs[j] = sqrt(xv[j]);
                x(_, j) = x(_, j) / xs[j];
                xv[j] = 1.0;
            } else {
                xs[j] = 1.0;
            }
        }
    } else {
        for (int j = 0; j < nvar; j++) {
            NumericVector xj = v * x(_, j);
            xv[j] = std::inner_product(xj.begin(), xj.end(), xj.begin(), 0.0);
            if (isd) {
                double xm = std::inner_product(v.begin(), v.end(), xj.begin(), 0.0);
                double xm_sqd = xm*xm;
                double vc = xv[j] - xm_sqd;
                xs[j] = sqrt(vc);
                x(_, j) = x(_, j) / xs[j];
                xv[j] = 1.0 + xm_sqd / vc;
            } else {
                xs[j] = 1.0;
            }
        }
    }
}

// [[Rcpp::export]]
void standardize_vec(NumericVector & y,
                     NumericVector & w,
                     double & ym,
                     double & ys,
                     bool & intr) {
    NumericVector v = sqrt(w);

    if (intr == true) {
        ym = std::inner_product(y.begin(), y.end(), w.begin(), 0.0);
        y = y - ym;
        ys = sqrt(std::inner_product(y.begin(), y.end(), y.begin(), 0.0) / y.size());
        y = y / ys;
    } else {
        NumericVector y_temp = v * y;
        double y_sqd = std::inner_product(y.begin(), y.end(), y.begin(), 0.0);
        double ym = std::inner_product(y.begin(), y.end(), v.begin(), 0.0);
        ys = sqrt(y_sqd - ym*ym);
        y = y / ys;
        ym = 0.0;
    }
}

//[[Rcpp:export]]
NumericVector compute_penalty(NumericVector & ulam,
                              int & nlam,
                              double & ptype,
                              double & pratio,
                              NumericVector & g,
                              NumericVector & cmult,
                              int & start,
                              int & stop) {

    NumericVector lambdas;
    if (ulam[0] == 0.0) {
        lambdas = NumericVector(nlam);
        lambdas[0] = 9.9e35;
        double max_pen = 0.0;
        for (int j = start; j < stop; j++) {
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

// [[Rcpp::export]]
void coord_desc(NumericMatrix & x,
                NumericVector & resid,
                NumericVector & ptype,
                int & no,
                int & nvar,
                int & nvar_total,
                NumericVector & cmult,
                NumericVector & upper_cl,
                NumericVector & lower_cl,
                int & ne,
                int & nx,
                NumericVector cur_lam,
                double & thr,
                int & maxit,
                NumericVector & xv,
                NumericMatrix & coef,
                NumericVector & b,
                NumericVector & g,
                NumericVector & rsq,
                double & rsq_cur,
                IntegerVector & mm,
                double & errcode,
                int & nlp,
                int & idx_lam) {

    int nin = 0;
    bool jerr = 0;
    NumericVector lasso_part(nvar_total);
    NumericVector ridge_part(nvar_total);

    for (int k = 0; k < nvar; k++) {
        lasso_part[k] = cmult[k] * ptype[k] * cur_lam[0];
        ridge_part[k] = xv[k] + cmult[k] * (1 - ptype[k]) * cur_lam[0];
    }

    for (int k = nvar; k < nvar_total; k++) {
            lasso_part[k] = cmult[k] * ptype[k] * cur_lam[1];
            ridge_part[k] = xv[k] + cmult[k] * (1 - ptype[k]) * cur_lam[1];
    }

    bool kkt_satisfied = 0, converge = 0;
    while(!kkt_satisfied & !jerr) {
        while(!converge & !jerr) {
            double dlx = 0.0;
            for (int k = 0; k < nvar_total; k++) {
                double gk = 0.0;
                for (int j = 0; j < no; j++) {
                    gk += (x(j, k) * resid[j]) / no;
                }
                double bk = b[k];
                double u = gk + bk * xv[k];
                double v = std::abs(u) - lasso_part[k];
                if (v > 0.0) {
                    b[k] = std::max(lower_cl[k], std::min(upper_cl[k], copysign(v, u) / ridge_part[k]));
                } else {
                    b[k] = 0.0;
                }
                if (b[k] != bk) {
                    if (mm[k] == 0) {
                        nin += 1;
                        mm[k] = nin;
                    }
                    if (nin <= nx) {
                        double del = b[k] - bk;
                        for (int j = 0; j < no; j++) {
                            resid[j] -= del * x(j, k);
                        }
                        rsq_cur += del * (2.0 * gk - del * xv[k]);
                        dlx = std::max(xv[k] * del * del, dlx);
                    } else {
                        jerr = 1;
                        errcode = -10000 - k;
                    }
                }
            }
            nlp += 1;
            if (nlp > maxit) {
                jerr = 1;
                errcode = -10000;
            }
            if (dlx < thr) {
                converge = 1;
            }
        }
        kkt_satisfied = 1;
    }
    coef(_, idx_lam) = b;
    rsq[idx_lam] = rsq_cur;
}

void unstandardize_coef(NumericMatrix & coef,
                        arma::mat & ext_,
                        int & nvar,
                        int & nvar_total,
                        int & nlam_total,
                        arma::vec & xm,
                        arma::vec & xs,
                        double & ys,
                        NumericVector & a0) {

    for (int i = 0; i < nvar; i++) {
        for (int j = 0; j < nlam_total; j++) {
            double z_alpha = a0[j];
            for (int k = nvar + 1; k < nvar_total; k++) {
                z_alpha += (ext_(i, k - nvar - 1) - xm[k]) / xs[k] * coef(k, j);
            }
            coef(i, j) = ys * (z_alpha + coef(i, j)) / xs[i];
        }
    }

}

// [[Rcpp::export]]
List gaussian_fit(int ka,
                  NumericVector ptype,
                  int nobs,
                  int nvar,
                  int nvar_ext,
                  arma::mat x_,
                  NumericVector y_,
                  arma::mat ext_,
                  NumericVector w,
                  NumericVector cmult,
                  NumericVector lower_cl,
                  NumericVector upper_cl,
                  int ne,
                  int ne_ext,
                  int nx,
                  int nx_ext,
                  int nlam,
                  int nlam_ext,
                  double pratio,
                  double pratio_ext,
                  NumericVector ulam_,
                  NumericVector ulam_ext_,
                  double thr,
                  int maxit,
                  bool isd,
                  bool isd_ext,
                  bool intr,
                  bool intr_ext) {

    // initialize internal variables
    int nvar_total;
    int nx_total;
    int ne_total;

    if (intr_ext == true) {
        nvar_total = nvar + nvar_ext + 1;
        nx_total = nx + nx_ext + 1;
        ne_total = ne + ne_ext + 1;
    } else {
        nvar_total = nvar + nvar_ext;
        nx_total = nx + nx_ext;
        ne_total = ne + ne_ext;
    }

    double ym;
    double ys;
    int nlam_total = nlam * nlam_ext;
    arma::vec xm(nvar_total);
    arma::vec xv(nvar_total);
    arma::vec xs(nvar_total);

    // copies of R objects to avoid modifying original
    NumericVector y = Rcpp::clone(y_);
    NumericVector wgt = w / sum(w);
    NumericVector ulam = Rcpp::clone(ulam_);
    NumericVector ulam_ext = Rcpp::clone(ulam_ext_);
    arma::vec w_arma = as<arma::vec>(wgt);

    // determine non-constant variables

    // Create single level regression matrix
    NumericMatrix xnew = create_data(nobs, nvar, nvar_total, x_, ext_, w_arma, isd, isd_ext, intr, intr_ext, xm, xv, xs);
    NumericVector xvnew = wrap(xv);

    // determine which variables are non-constant

    // standardize y, confidence limits --> user penalties??
    standardize_vec(y, wgt, ym, ys, intr);
    lower_cl = lower_cl / ys;
    upper_cl = upper_cl / ys;
    //ulam = ulam_ / ys;
    //ulam_ext = ulam_ext_ / ys;

    if (isd == true) {
        for (int i = 0; i < nvar; i++) {
            if (lower_cl[i] != R_NegInf) {
                lower_cl[i] = lower_cl[i] * xs[i];
            }
            if (upper_cl[i] != R_PosInf) {
                upper_cl[i] = upper_cl[i] * xs[i];
            }
        }
    }
    if (isd_ext == true) {
        for (int i = nvar; i < nvar_total; i++) {
            if (lower_cl[i] != R_NegInf) {
                lower_cl[i] = lower_cl[i] * xs[i];
            }
            if (upper_cl[i] != R_PosInf) {
                upper_cl[i] = upper_cl[i] * xs[i];
            }
        }
    }

    // initialize objects to hold fitting results
    NumericVector coef0(nlam_total, 0.0);
    NumericMatrix coef(nx_total, nlam_total);
    NumericVector rsq(nlam_total);

    // ---------run coordinate descent for all penalties ----------
    NumericVector g(nvar_total, 0.0);
    IntegerVector mm(nvar_total, 0);

    // compute penalty paths
    for (int k = 0; k < nvar_total; k++) {
        g[k] = std::abs(std::inner_product(xnew(_, k).begin(), xnew(_, k).end(), y.begin(), 0.0)) / nobs;
    }

    int start = 0;
    NumericVector lam_path = compute_penalty(ulam, nlam, ptype[1], pratio, g, cmult, start, nvar);
    start = nvar + 1;
    NumericVector lam_path_ext = compute_penalty(ulam_ext, nlam_ext, ptype[nvar + 1], pratio_ext, g, cmult, start, nvar_total);

    // loop through all penalty combinations
    int nlp = 0;
    bool flag = false;
    int idx_lam = 0;
    double errcode = 0.0;
    NumericVector lam_cur(2, 0.0);

    double rsq_outer = 0.0;
    double rsq_inner = 0.0;

    NumericVector outer_resid = Rcpp::clone(y);
    NumericVector inner_resid(nobs);

    NumericVector coef_outer(nvar_total, 0.0);
    NumericVector coef_inner(nvar_total, 0.0);

    for (int m = 0; m < nlam; m++) {

        lam_cur[0] = lam_path[m];

        for (int m2 = 0; m2 < nlam_ext; m2++) {

            lam_cur[1] = lam_path_ext[m2];

            if (m2 == 0) {
                coord_desc(xnew, outer_resid, ptype, nobs, nvar, nvar_total, cmult, upper_cl, lower_cl, ne_total, nx_total,
                           lam_cur, thr, maxit, xvnew, coef, coef_outer, g, rsq, rsq_outer, mm, errcode, nlp, idx_lam);

                inner_resid = Rcpp::clone(outer_resid);
                coef_inner = Rcpp::clone(coef_outer);
                rsq_inner = rsq_outer;

            } else {
                coord_desc(xnew, inner_resid, ptype, nobs, nvar, nvar_total, cmult, upper_cl, lower_cl, ne_total, nx_total,
                           lam_cur, thr, maxit, xvnew, coef, coef_inner, g, rsq, rsq_inner, mm, errcode, nlp, idx_lam);
            }

            // stop if max iterations or no change in r-squared
            //if ((rsq_inner - rsq_old) < 1e-05 * rsq_inner || rsq_inner > 0.999 || errcode != 0.0) {
            //    flag = true;
            //} else {
            //    rsq_old = rsq_inner;
            //    ++idx_lam;
            //}

            ++idx_lam;
        }
        if (flag) {
            break;
        }
    }

    // unstandardize predictor variables
    NumericVector b0(nlam_total, 0.0);
    NumericVector a0(nlam_total, 0.0);

    if (intr_ext == true) {
        a0 = coef(nvar, _);
    }
    unstandardize_coef(coef, ext_, nvar, nvar_total, nlam_total, xm, xs, ys, a0);

    // unstandardize external variables
    for (int i = (nvar + 1); i < nvar_total; i++) {
        for(int j = 0; j < nlam_total; j++) {
                coef(i, j) = ys * coef(i, j) / xs[i];
        }
    }
    a0 = ys * a0;

    // return model fit for all penalty combinations
    return Rcpp::List::create(Named("lam") = lam_path,
                              Named("lam_ext") = lam_path_ext,
                              Named("coef") = coef,
                              Named("nlp") = nlp,
                              Named("pratioext") = pratio_ext,
                              Named("xnew") = xnew,
                              Named("rsq") = rsq);
}
