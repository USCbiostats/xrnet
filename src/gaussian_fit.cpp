#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericMatrix create_data(const int & nobs,
                          const int & nvar,
                          const int & nvar_total,
                          const arma::mat & x,
                          const arma::mat & ext,
                          const arma::vec & w,
                          const bool & isd,
                          const bool & isd_ext,
                          const bool & intr,
                          const bool & intr_ext,
                          arma::vec & xm,
                          arma::vec & xv,
                          arma::vec & xs) {

    // initialize new design matrix
    arma::mat xnew(nobs, nvar_total);

    // standardize x (predictor variables)
    if (intr == true) {
        if (isd == true) {
            for (int j = 0; j < nvar; j++) {
                xm[j] = arma::dot(w, x.col(j));
                xs[j] = sqrt(arma::dot(x.col(j) - xm[j], (x.col(j) - xm[j]) % w));
                xnew.col(j) = (x.col(j) - xm[j]) / xs[j];
            }
        }
        else {
            for (int j = 0; j < nvar; j++) {
                xm[j] = arma::dot(x.col(j), w);
                xnew.col(j) = x.col(j) - xm[j];
                xv[j] = arma::dot(xnew.col(j), xnew.col(j) % w);
            }
        }
    }
    else {
        if (isd == true) {
            for (int j = 0; j < nvar; j++) {
                double xm_j = arma::dot(w, x.col(j));
                double vc = arma::dot(x.col(j) - xm_j, (x.col(j) - xm_j) % w);
                xs[j] = sqrt(vc);
                xnew.col(j) = x.col(j) / xs[j];
                xv[j] = 1.0 + xm_j * xm_j / vc;
            }
        }
        else {
            for (int j = 0; j < nvar; j++) {
                xnew.col(j) = x.col(j);
                xv[j] = arma::dot(x.col(j), x.col(j) % w);
            }
        }
    }

    // standarize ext (external data variables)
    const arma::mat xsub = xnew.submat(0, 0, nobs - 1, nvar - 1);

    if (intr_ext == true) {

        // create intercept column for external data
        xnew.col(nvar) = xsub * arma::ones<arma::mat>(nvar, 1);
        xv[nvar] = arma::dot(xnew.col(nvar), xnew.col(nvar)) / nobs;
        xs[nvar] = 1.0;

        if (isd_ext == true) {
            for (int j = nvar + 1; j < nvar_total; j++) {
                xm[j] = arma::mean(ext.col(j - nvar - 1));
                xs[j] = sqrt(arma::dot(ext.col(j - nvar - 1) - xm[j], ext.col(j - nvar - 1) - xm[j]) / nvar);
                xnew.col(j) = xsub * ((ext.col(j - nvar - 1) - xm[j]) / xs[j]);
                xv[j] = arma::dot(xnew.col(j), xnew.col(j)) / nobs;
            }
        }
        else {
            for (int j = nvar + 1; j < nvar_total; j++) {
                xm[j] = arma::mean(ext.col(j - nvar - 1));
                xnew.col(j) = xsub * (ext.col(j - nvar - 1)  - xm[j]);
                xv[j] = arma::dot(xnew.col(j), xnew.col(j)) / nobs;
            }
        }
    }
    else {
        if (isd_ext == true) {
            for (int j = nvar; j < nvar_total; j++) {
                double xm_j = arma::mean(ext.col(j - nvar));
                xs[j] = sqrt(arma::dot(ext.col(j - nvar) - xm_j, ext.col(j - nvar) - xm_j) / nvar);
                xnew.col(j) = xsub * (ext.col(j - nvar) / xs[j]);
                xv[j] = arma::dot(xnew.col(j), xnew.col(j)) / nobs;
            }
        }
        else {
            for (int j = nvar; j < nvar_total; j++) {
                xnew.col(j) = xsub * ext.col(j - nvar);
                xv[j] = arma::dot(xnew.col(j), xnew.col(j)) / nobs;
            }
        }
    }
    return(wrap(xnew));
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
        double y_sqd = std::inner_product(y.begin(), y.end(), y.begin(), 0.0) / y.size();
        double ym = std::inner_product(y.begin(), y.end(), w.begin(), 0.0);
        ys = sqrt(y_sqd - ym*ym);
        y = y / ys;
        ym = 0.0;
    }
}

//[[Rcpp:export]]
NumericVector compute_penalty(NumericVector & ulam,
                              const int & nlam,
                              const double & ptype,
                              const double & pratio,
                              NumericVector & g,
                              NumericVector & cmult,
                              const int & start,
                              const int & stop) {

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
void coord_desc(const NumericMatrix & x,
                NumericVector & resid,
                const NumericVector & ptype,
                const int & no,
                const int & nvar,
                const int & nvar_total,
                const NumericVector & cmult,
                const NumericVector & upper_cl,
                const NumericVector & lower_cl,
                const int & ne,
                const int & nx,
                NumericVector cur_lam,
                const double & thr,
                const int & maxit,
                const NumericVector & xv,
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
                        const arma::mat & ext_,
                        const int & nvar,
                        const int & nvar_total,
                        int & nlam_total,
                        const arma::vec & xm,
                        const arma::vec & xs,
                        double & ys,
                        NumericVector & a0,
                        bool & intr_ext) {

    int start_idx = nvar;
    if (intr_ext == true) {
        start_idx += 1;
    }

    for (int i = 0; i < nvar; i++) {
        for (int j = 0; j < nlam_total; j++) {
            double z_alpha = a0[j];
            for (int k = start_idx; k < nvar_total; k++) {
                z_alpha += (ext_(i, k - start_idx) - xm[k]) / xs[k] * coef(k, j);
            }
            coef(i, j) = ys * (z_alpha + coef(i, j)) / xs[i];
        }
    }

}

// [[Rcpp::export]]
List gaussian_fit(const int ka,
                  NumericVector ptype,
                  const int nobs,
                  const int nvar,
                  const int nvar_ext,
                  const arma::mat & x_,
                  NumericVector y_,
                  const arma::mat & ext_,
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
    arma::vec xm(nvar_total, arma::fill::zeros);
    arma::vec xv(nvar_total, arma::fill::ones);
    arma::vec xs(nvar_total, arma::fill::ones);

    // copies of R objects to avoid modifying original
    NumericVector y = Rcpp::clone(y_);
    NumericVector wgt = w / sum(w);
    NumericVector ulam = Rcpp::clone(ulam_);
    NumericVector ulam_ext = Rcpp::clone(ulam_ext_);
    const arma::vec w_arma = as<arma::vec>(wgt);

    // determine non-constant variables

    // Create single level regression matrix
    NumericMatrix xnew = create_data(nobs, nvar, nvar_total, x_, ext_, w_arma, isd, isd_ext, intr, intr_ext, xm, xv, xs);
    NumericVector xvnew = wrap(xv);

    // determine which variables are non-constant

    // standardize y, confidence limits, user penalties
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
    //bool flag = false;
    int idx_lam = 0;
    double errcode = 0.0;
    NumericVector lam_cur(2, 0.0);

    double rsq_outer = 0.0;
    double rsq_inner = 0.0;
    double rsq_old = 0.0;

    NumericVector outer_resid = Rcpp::clone(y);
    NumericVector inner_resid(nobs);

    NumericVector coef_outer(nvar_total, 0.0);
    NumericVector coef_inner(nvar_total, 0.0);

    for (int m = 0; m < nlam; m++) {
        lam_cur[0] = lam_path[m];

        for (int m2 = 0; m2 < nlam_ext; m2++) {
            lam_cur[1] = lam_path_ext[m2];

            if (m2 == 0) {
                coord_desc(xnew, outer_resid, ptype, nobs, nvar, nvar_total,
                           cmult, upper_cl, lower_cl, ne_total, nx_total,
                           lam_cur, thr, maxit, xvnew, coef, coef_outer, g,
                           rsq, rsq_outer, mm, errcode, nlp, idx_lam);

                inner_resid = Rcpp::clone(outer_resid);
                coef_inner = Rcpp::clone(coef_outer);
                rsq_inner = rsq_outer;
                rsq_old = rsq_outer;

            }
            else {
                coord_desc(xnew, inner_resid, ptype, nobs, nvar, nvar_total,
                           cmult, upper_cl, lower_cl, ne_total, nx_total,
                           lam_cur, thr, maxit, xvnew, coef, coef_inner, g,
                           rsq, rsq_inner, mm, errcode, nlp, idx_lam);

                //stop if max iterations or no change in r-squared
                if ((rsq_inner - rsq_old) < (1e-05 * rsq_inner) || rsq_inner > 0.999 || errcode != 0.0) {
                    break;
                }
                else {
                    rsq_old = rsq_inner;
                    ++idx_lam;
                }
            }
        }
    }

    // unstandardize predictor variables
    NumericVector b0(nlam_total, 0.0);
    NumericVector a0(nlam_total, 0.0);

    if (intr_ext == true) {
        a0 = coef(nvar, _);
    }
    unstandardize_coef(coef, ext_, nvar, nvar_total, nlam_total, xm, xs, ys, a0, intr_ext);

    // unstandardize external variables
    for (int i = nvar; i < nvar_total; i++) {
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
                              Named("rsq") = rsq);
}
