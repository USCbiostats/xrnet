#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/*
// [[Rcpp::export]]
NumericMatrix create_data(int & nobs,
                               int & nvar,
                               int & nvar_total,
                               NumericMatrix & x,
                               NumericMatrix & ext,
                               NumericVector & w,
                               bool & isd,
                               bool & isd_ext,
                               bool & intr,
                               bool & intr_ext,
                               NumericVector & xm,
                               NumericVector & xv,
                               NumericVector & xs) {

    NumericVector v = sqrt(w);
    NumericMatrix xnew(nobs, nvar_total);

    if (intr == true) {
        for (int j = 0; j < nvar; j++) {
            xm[j] = std::inner_product(x(_, j).begin(), x(_, j).end(), w.begin(), 0.0);
            xnew(_, j) = x(_, j) - xm[j];
            xv[j] = std::inner_product(x(_, j).begin(), x(_, j).end(), x(_, j).begin(), 0.0) / nobs;
            if (isd == true) {
                xs[j] = sqrt(xv[j]);
                xnew(_, j) = xnew(_, j) / xs[j];
                xv[j] = 1.0;
            } else {
                xs[j] = 1.0;
            }
        }
    } else {
        for (int j = 0; j < nvar; j++) {
            NumericVector xj = v * x(_, j);
            xv[j] = std::inner_product(xj.begin(), xj.end(), xj.begin(), 0.0);
            if (isd == true) {
                double xm = std::inner_product(v.begin(), v.end(), xj.begin(), 0.0);
                double vc = xv[j] - xm*xm;
                xs[j] = sqrt(vc);
                xnew(_, j) = x(_, j) / xs[j];
                xv[j] = 1.0 + xm*xm / vc;
            } else {
                xs[j] = 1.0;
            }
        }
    }

    if (intr_ext == true) {

        for (int j = nvar; j < nvar_total; j++) {
            if (j == nvar) {

            }
        }

    } else {

    }

    return(xnew);
}

*/

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
        for (int j = 0; j < nvar; j++) {
            xnew.col(j) = x.col(j) - xm[j];
            xv[j] = arma::dot(xnew.col(j), xnew.col(j)) / nobs;
            if (isd == true) {
                xs[j] = sqrt(xv[j]);
                xnew.col(j) = xnew.col(j) / xs[j];
                xv[j] = 1.0;
            } else {
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

    arma::vec ext_temp(nvar);

    if (intr_ext) {
        xnew.col(nvar) = xnew.submat(0, 0, nobs - 1, nvar - 1) * arma::ones<arma::mat>(nvar, 1);
        xm.subvec(nvar + 1, nvar_total - 1) = arma::mean(ext).t();
        for (int j = nvar + 1; j < nvar_total; j++) {
            ext_temp = ext.col(j - nvar - 1) - xm[j];
            xv[j] = arma::dot(ext_temp, ext_temp) / nvar;
            if (isd_ext == true) {
                xs[j] = sqrt(xv[j]);
                ext_temp = ext_temp / xs[j];
                xv[j] = 1.0;
            } else {
                xs[j] = 1.0;

            }
            xnew.col(j) = xnew.submat(0, 0, nobs - 1, nvar - 1) * ext_temp;
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
        lambdas[1] = alf * max_pen;
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
                int & nv,
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
    NumericVector lasso_part(nv);
    NumericVector ridge_part(nv);

    for (int k = 0; k < nv; k++) {
        if (ptype[k] <= 1) {
            lasso_part[k] = cmult[k] * ptype[k] * cur_lam[k];
            ridge_part[k] = xv[k] + cmult[k] * (1 - ptype[k]) * cur_lam[k];
        }
    }

    bool kkt_satisfied = 0, converge = 0;
    while(!kkt_satisfied & !jerr) {
        while(!converge & !jerr) {
            double dlx = 0.0;
            for (int k = 0; k < nv; k++) {
                double gk = std::inner_product(x(_, k).begin(), x(_, k).end(), resid.begin(), 0.0) / no;
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
                        rsq_cur = rsq_cur + del * (2.0 * gk - del * xv[k]);
                        resid = resid - del * x(_, k);
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

// [[Rcpp::export]]
List gaussian_fit(int ka,
                  NumericVector ptype,
                  int nobs,
                  int nvar,
                  int nvar_ext,
                  NumericMatrix x_,
                  NumericVector y_,
                  NumericMatrix ext_,
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
    NumericVector x_mean(nvar, 0.0);
    NumericVector ext_mean(nvar_ext, 0.0);
    NumericVector x_var(nvar, 0.0);
    NumericVector ext_var(nvar_ext, 0.0);
    NumericVector x_sd(nvar, 0.0);
    NumericVector ext_sd(nvar_ext, 0.0);

    // copies of R objects to avoid modifying original
    NumericMatrix x = Rcpp::clone(x_);
    NumericVector y = Rcpp::clone(y_);
    NumericVector wgt = w / sum(w);
    NumericVector w_ext(nvar, 1.0 / nvar);
    NumericMatrix ext = Rcpp::clone(ext_);
    NumericVector ulam = Rcpp::clone(ulam_);
    NumericVector ulam_ext = Rcpp::clone(ulam_ext_);

    // determine which variables are non-constant

    // standarize x, ext, y, confidence limits --> standardize user penalties??
    standardize_mat(nobs, nvar, x, wgt, isd, intr, x_mean, x_var, x_sd);
    standardize_mat(nvar, nvar_ext, ext, w_ext, isd_ext, intr_ext, ext_mean, ext_var, ext_sd);
    standardize_vec(y, wgt, ym, ys, intr);
    lower_cl = lower_cl / ys;
    upper_cl = upper_cl / ys;

    if (isd == true) {
        for (int i = 0; i < nvar; i++) {
            if (lower_cl[i] != R_NegInf) {
                lower_cl[i] = lower_cl[i] * x_sd[i];
            }
            if (upper_cl[i] != R_PosInf) {
                upper_cl[i] = upper_cl[i] * x_sd[i];
            }
        }
    }
    if (isd_ext == true) {
        for (int i = nvar; i < nvar_total; i++) {
            if (lower_cl[i] != R_NegInf) {
                lower_cl[i] = lower_cl[i] * ext_sd[i];
            }
            if (upper_cl[i] != R_PosInf) {
                upper_cl[i] = upper_cl[i] * ext_sd[i];
            }
        }
    }

    // create modified design matrix and parameters for input to coordinate descent
    NumericMatrix xnew;
    if (intr_ext == true) {
        xnew = wrap(arma::join_rows(as<arma::mat>(x), as<arma::mat>(x) * join_rows(arma::ones<arma::mat>(nvar, 1), as<arma::mat>(ext))));
    } else {
        xnew = wrap(arma::join_rows(as<arma::mat>(x), as<arma::mat>(x) * as<arma::mat>(ext)));
    }

    NumericVector xv(nvar_total, 0.0);
    std::copy(x_var.begin(), x_var.end(), xv.begin());
    for (int k = nvar; k < nvar_total; k++) {
        double xm = std::inner_product(xnew(_, k).begin(), xnew(_, k).end(), wgt.begin(), 0.0);
        double xm2 = std::inner_product(xnew(_, k).begin(), xnew(_, k).end(), xnew(_, k).begin(), 0.0) / nobs;
        xv[k] = xm2 - xm*xm;
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
    NumericVector lam_path = compute_penalty(ulam, nlam, pratio, g, cmult, start, nvar);
    NumericVector lam_path_ext = compute_penalty(ulam_ext, nlam_ext, pratio_ext, g, cmult, nvar, nvar_total);

    // loop through all penalty combinations
    int nlp = 0;
    bool flag = false;
    int idx_lam = 0;
    double errcode = 0.0;
    NumericVector lam_cur(nvar_total, 0.0);

    double rsq_outer = 0.0;
    double rsq_inner = 0.0;
    //double rsq_old = 0.0;

    NumericVector outer_resid = Rcpp::clone(y);
    NumericVector inner_resid(nobs);

    NumericVector coef_outer(nvar_total, 0.0);
    NumericVector coef_inner(nvar_total, 0.0);

    for (int m = 0; m < nlam; m++) {
        for (int m2 = 0; m2 < nlam_ext; m2++) {

            for (int j = 0; j < nvar; j++) {
                lam_cur[j] = lam_path[m];
            }

            for (int j = nvar; j < nvar_total; j++) {
                lam_cur[j] = lam_path_ext[m2];
            }

            if (m2 == 0) {
                coord_desc(xnew, outer_resid, ptype, nobs, nvar_total, cmult, upper_cl, lower_cl, ne_total, nx_total,
                           lam_cur, thr, maxit, xv, coef, coef_outer, g, rsq, rsq_outer, mm, errcode, nlp, idx_lam);

                inner_resid = outer_resid;
                coef_inner = coef_outer;
                rsq_inner = rsq_outer;

            } else {
                coord_desc(xnew, inner_resid, ptype, nobs, nvar_total, cmult, upper_cl, lower_cl, ne_total, nx_total,
                           lam_cur, thr, maxit, xv, coef, coef_inner, g, rsq, rsq_inner, mm, errcode, nlp, idx_lam);
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

    // return model fit for all penalty combinations
    return Rcpp::List::create(Named("lam") = lam_path,
                              Named("lam_ext") = lam_path_ext,
                              Named("coef") = coef,
                              Named("nlp") = nlp,
                              Named("pratioext") = pratio_ext);
}
