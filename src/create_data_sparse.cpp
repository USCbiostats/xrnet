#include <RcppArmadillo.h>
#include <cmath>
#include <numeric>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/*
* Generates matrix [x | unpen | x*external]
* external data is sparse matrix
*/

// [[Rcpp::export]]
double mean_sparse(const arma::sp_mat & x, int & col_cur, const int & nvar) {
    double xsum = 0.0;
    for(arma::sp_mat::const_iterator it = x.begin_col(col_cur); it != x.end_col(col_cur); ++it) {
        xsum += *it;
    }
    return xsum / nvar;
}

// [[Rcpp::export]]
double sd_sparse(const arma::sp_mat & x, int & col_cur, double & xm, const int & nvar) {
    double xsum2 = 0.0;
    int nnzero = 0;
    for(arma::sp_mat::const_iterator it = x.begin_col(col_cur); it != x.end_col(col_cur); ++it) {
        xsum2 += (*it - xm) * (*it - xm);
        nnzero++;
    }
    xsum2 += (nvar - nnzero) * xm * xm;
    return std::sqrt(xsum2 / nvar);
}

// [[Rcpp::export]]
arma::mat create_data_sparse(const int & nobs,
                             const int & nvar,
                             const int & nvar_ext,
                             const int & nvar_unpen,
                             const int & nvar_total,
                             const arma::mat & x,
                             const arma::sp_mat & ext,
                             const arma::mat & unpen,
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
    arma::mat xsub(xnew.memptr(), nobs, nvar, false, false);

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
    if (nvar_ext > 0) {
        if (intr_ext) {
            if (isd_ext) {
                for (int j = 0; j < nvar_ext; ++j) {
                    xm[xnew_col] = mean_sparse(ext, j, nvar);
                    xs[xnew_col] = sd_sparse(ext, j, xm[xnew_col], nvar);
                    xnew.col(xnew_col) = xsub * ((arma::vec(ext.col(j)) - xm[xnew_col]) / xs[xnew_col]);
                    xv[xnew_col] = arma::dot(xnew.unsafe_col(xnew_col), xnew.unsafe_col(xnew_col)) / nobs;
                    ++xnew_col;
                }
            }
            else {
                for (int j = 0; j < nvar_ext; ++j) {
                    xm[xnew_col] = mean_sparse(ext, j, nvar);
                    xnew.col(xnew_col) = xsub * (arma::vec(ext.col(j)) - xm[xnew_col]);
                    xv[xnew_col] = arma::dot(xnew.unsafe_col(xnew_col), xnew.unsafe_col(xnew_col)) / nobs;
                    ++xnew_col;
                }
            }
        }
        else {
            if (isd_ext) {
                for (int j = 0; j < nvar_ext; ++j) {
                    double xm_j = mean_sparse(ext, j, nvar);
                    xs[xnew_col] = sd_sparse(ext, j, xm_j, nvar);
                    xnew.col(xnew_col) = xsub * (arma::vec(ext.col(j)) / xs[xnew_col]);
                    xv[xnew_col] = arma::dot(xnew.unsafe_col(xnew_col), xnew.unsafe_col(xnew_col)) / nobs;
                    ++xnew_col;
                }
            }
            else {
                for (int j = 0; j < nvar_ext; ++j) {
                    xnew.col(xnew_col) = xsub * arma::vec(ext.col(j));
                    xv[xnew_col] = arma::dot(xnew.unsafe_col(xnew_col), xnew.unsafe_col(xnew_col)) / nobs;
                    ++xnew_col;
                }
            }
        }
    }
    return(xnew);
}
