#ifndef HIERR_UTILS_H
#define HIERR_UTILS_H

#include <RcppArmadillo.h>

template <class T>
inline T sgn(T & v) {
    return (v > T(0)) - (v < T(0));
}

template <class T>
void save_results(arma::vec coef_cur,
                  arma::vec & a0_coef,
                  arma::mat & alphas_coef,
                  arma::vec & b0_coef,
                  arma::mat & betas_coef,
                  arma::mat & gammas_coef,
                  const T & ext,
                  const int & nvar,
                  const int & nvar_unpen,
                  const int & nv_x,
                  const int & ext_start,
                  const int & nvar_ext,
                  const int & nvar_total,
                  const arma::vec & xm,
                  const arma::vec & xs,
                  const double & ym,
                  const double & ys,
                  const bool & intr,
                  const bool & intr_ext,
                  const int & idx_lam) {

    // unstandardize variables
    coef_cur = ys * (coef_cur / xs);

    // get external coef
    if (nvar_ext > 0) {
        alphas_coef.col(idx_lam) = coef_cur.tail(nvar_ext);
    }

    // unstandardize predictors w/ external data
    if (nvar_ext + intr_ext > 0) {
        arma::vec z_alpha(nvar);
        if (intr_ext) {
            z_alpha.fill(coef_cur[nv_x]);
        } else {
            z_alpha.fill(0.0);
        }
        if (nvar_ext > 0) {
            z_alpha += arma::vec(ext * coef_cur.tail(nvar_ext));
        }
        betas_coef.col(idx_lam) = z_alpha / xs.head(nvar) + coef_cur.head(nvar);
    } else {
        betas_coef.col(idx_lam) = coef_cur.head(nvar);
    }

    // unstandardize predictors w/o external data
    if (nvar_unpen > 0) {
        gammas_coef.col(idx_lam) = coef_cur.subvec(nvar, nv_x - 1);
    }

    // compute 2nd level intercepts
    if (intr_ext) {
        if (nvar_ext > 0) {
            a0_coef[idx_lam] = arma::mean(betas_coef.col(idx_lam)) - arma::dot(xm.tail(nvar_ext), alphas_coef.col(idx_lam));
        } else {
            a0_coef[idx_lam] = arma::mean(betas_coef.col(idx_lam));
        }
    }

    // compute 1st level intercepts
    if (intr) {
        b0_coef[idx_lam] = ym - arma::dot(xm.head(nvar), betas_coef.unsafe_col(idx_lam));
        if (nvar_unpen > 0) {
            b0_coef[idx_lam] -= arma::dot(xm.subvec(nvar, nv_x - 1), gammas_coef.unsafe_col(idx_lam));
        }
    }
}

int countNonzero(const arma::vec & x,
                 const int & start,
                 const int & end);

void updatePenalty(arma::vec & l1,
                   arma::vec & l2,
                   const arma::vec & pind,
                   const arma::vec & cmult,
                   const arma::vec & xv,
                   const double & lam_new,
                   const int & start,
                   const int & stop);

void updateStrong(Rcpp::LogicalVector & strong,
                  const arma::vec & g,
                  const arma::vec & ptype,
                  const arma::vec & cmult,
                  const Rcpp::NumericVector & lam_cur,
                  const Rcpp::NumericVector & lam_prev,
                  const Rcpp::NumericVector & qnt,
                  const Rcpp::IntegerVector & blkend);

void compute_penalty(Rcpp::NumericVector & lambdas,
                     Rcpp::NumericVector & ulam,
                     const double & ptype,
                     const double & pratio,
                     const arma::vec & g,
                     const arma::vec & cmult,
                     const int & start,
                     const int & stop);

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
                  const int & ext_start);

void compute_coef_sparse(arma::mat & coef,
                         const arma::sp_mat & ext_,
                         const int & nvar,
                         const int & nvar_ext,
                         const int & nvar_total,
                         const int & nlam_total,
                         const arma::vec & xm,
                         const arma::vec & xs,
                         const double & ys,
                         const arma::vec & a0,
                         const bool & intr_ext,
                         const int & ext_start);

void standardize_vec(arma::vec & y,
                     const arma::vec & w,
                     double & ym,
                     double & ys,
                     const bool & intr);

#endif
