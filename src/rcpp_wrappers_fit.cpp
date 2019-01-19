// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <string.h>
#include "DataFunctions.h"
#include "Hierr.h"
#include "HierrUtils.h"
#include "CoordDescTypes.h"
#include "BinomialSolver.h"

/*
* Rcpp wrapper to fit model when x is big.matrix or
* filetracked.big.matrix and external is dense matrix
*/

// [[Rcpp::export]]
Rcpp::List fitDenseDense(SEXP xbig,
                         const Eigen::Map<Eigen::VectorXd> y,
                         const Eigen::Map<Eigen::MatrixXd> ext,
                         const Eigen::Map<Eigen::MatrixXd> fixed,
                         Eigen::VectorXd weights_user,
                         const Rcpp::LogicalVector & intr,
                         const Rcpp::LogicalVector & stnd,
                         const Eigen::Map<Eigen::VectorXd> penalty_type,
                         const Eigen::Map<Eigen::VectorXd> cmult,
                         const Eigen::Map<Eigen::VectorXd> quantiles,
                         const Rcpp::IntegerVector & num_penalty,
                         const Rcpp::NumericVector & penalty_ratio,
                         const Eigen::Map<Eigen::VectorXd> penalty_user,
                         const Eigen::Map<Eigen::VectorXd> penalty_user_ext,
                         Eigen::VectorXd lower_cl,
                         Eigen::VectorXd upper_cl,
                         const std::string & family,
                         const double & thresh,
                         const int & maxit,
                         const int & ne,
                         const int & nx) {

    // get pointer to big.matrix for X and map to eigen matrix
    Rcpp::XPtr<BigMatrix> xptr(xbig);
    Eigen::Map<const Eigen::MatrixXd> x((const double *)xptr->matrix(),
                                        xptr->nrow(),
                                        xptr->ncol());

    // initialize objects to hold means, variances, sds of all variables
    const int n = x.rows();
    const int nv_x = x.cols();
    const int nv_fixed = fixed.size() == 0 ? 0 : fixed.cols();
    const int nv_ext = ext.size() == 0 ? 0 : ext.cols();
    const int nv_total = nv_x + nv_fixed + intr[1] + nv_ext;
    Eigen::VectorXd xm = Eigen::VectorXd::Constant(nv_total, 0.0);
    Eigen::VectorXd xv = Eigen::VectorXd::Constant(nv_total, 1.0);
    Eigen::VectorXd xs = Eigen::VectorXd::Constant(nv_total, 1.0);

    // map to correct size of fixed matrix
    Eigen::Map<const Eigen::MatrixXd> fixedmap(fixed.data(), fixed.rows(), nv_fixed);

    // scale user weights
    weights_user.array() = weights_user.array() / weights_user.sum();

    // compute moments of matrices and create XZ (if external data present)
    compute_moments(x, weights_user, xm, xv, xs, true, stnd[0], 0);
    compute_moments(fixedmap, weights_user, xm, xv, xs, true, stnd[0], nv_x);
    const Eigen::MatrixXd xz = create_XZ(x, ext, xm, xv, xs, intr[1]);
    compute_moments(xz, weights_user, xm, xv, xs, false, stnd[1], nv_x + nv_fixed);

    // choose solver based on outcome
    std::unique_ptr<CoordSolver<MapMat> > solver;
    if (family == "gaussian") {
        solver.reset(new CoordSolver<MapMat>(y, x, fixedmap, xz, xm.data(), xv.data(), xs.data(),
                                             weights_user, intr[0], penalty_type.data(),
                                             cmult.data(), quantiles, upper_cl.data(),
                                             lower_cl.data(), ne, nx, thresh, maxit));

    }
    else if (family == "binomial") {
        solver.reset(new BinomialSolver<MapMat>(y, x, fixedmap, xz, xm.data(), xv.data(), xs.data(),
                                                weights_user, intr[0], penalty_type.data(),
                                                cmult.data(), quantiles, upper_cl.data(),
                                                lower_cl.data(), ne, nx, thresh, maxit));
    }

    // Object to hold results for all penalty combinations
    const int num_combn = num_penalty[0] * num_penalty[1];
    Hierr<MapMat, MapMat> estimates = Hierr<MapMat, MapMat>(n, nv_x, nv_fixed, nv_ext, nv_total,
                                                            intr[0], intr[1], ext.data(), xm.data(),
                                                            xs.data(), num_combn);

    // compute penalty path for 1st level variables
    Eigen::VectorXd path(num_penalty[0]);

    compute_penalty(path, penalty_user, penalty_type[0],
                    penalty_ratio[0], solver->getGradient(),
                    solver->getCmult(), 0, nv_x);

    // compute penalty path for 2nd level variables
    Eigen::VectorXd path_ext(num_penalty[1]);
    if (nv_ext > 0) {
        compute_penalty(path_ext, penalty_user_ext,
                        penalty_type[nv_x + nv_fixed + intr[1]],
                                    penalty_ratio[1], solver->getGradient(),
                                    solver->getCmult(), nv_x + nv_fixed + intr[1], nv_total);
    } else {
        path_ext[0] = 0.0;
    }

    // solve grid of penalties in decreasing order
    double b0_outer = solver->getBeta0();
    Eigen::VectorXd betas_outer = solver->getBetas();
    Eigen::VectorXd strong_sum = Eigen::VectorXd::Zero(num_combn);
    Eigen::VectorXd active_sum = Eigen::VectorXd::Zero(num_combn);

    int idx_pen = 0;
    for (int m = 0; m < num_penalty[0]; ++m) {
        solver->setPenalty(path[m], 0);
        for (int m2 = 0; m2 < num_penalty[1]; ++m2, ++idx_pen) {
            solver->setPenalty(path_ext[m2], 1);
            if (m2 == 0 && num_penalty[1] > 1) {
                solver->warm_start(y, b0_outer, betas_outer);
                solver->update_strong(path, path_ext, m, m2);
                solver->solve();
                b0_outer = solver->getBeta0();
                betas_outer = solver->getBetas();
            }
            else {
                solver->update_strong(path, path_ext, m, m2);
                solver->solve();
            }
            strong_sum[idx_pen] = sum(solver->getStrongSet());
            active_sum[idx_pen] = sum(solver->getActiveSet());
            estimates.add_results(solver->getBeta0(), solver->getBetas(), idx_pen);
        }
    }

    int status = 0;

    // collect results in list and return to R
    return Rcpp::List::create(Rcpp::Named("beta0") = estimates.getBeta0(),
                              Rcpp::Named("betas") = estimates.getBetas(),
                              Rcpp::Named("alpha0") = estimates.getAlpha0(),
                              Rcpp::Named("alphas") = estimates.getAlphas(),
                              Rcpp::Named("num_passes") = solver->getNumPasses(),
                              Rcpp::Named("penalty") = path,
                              Rcpp::Named("penalty_ext") = path_ext,
                              Rcpp::Named("strong_sum") = strong_sum,
                              Rcpp::Named("active_sum") = active_sum,
                              Rcpp::Named("status") = status);
}

/*
* Rcpp wrapper to fit model when x is sparse matrix
*  and external is dense matrix
*/

// [[Rcpp::export]]
Rcpp::List fitSparseDense(const Eigen::MappedSparseMatrix<double> x,
                          const Eigen::Map<Eigen::VectorXd> y,
                          const Eigen::Map<Eigen::MatrixXd> ext,
                          const Eigen::Map<Eigen::MatrixXd> fixed,
                          Eigen::VectorXd weights_user,
                          const Rcpp::LogicalVector & intr,
                          const Rcpp::LogicalVector & stnd,
                          const Eigen::Map<Eigen::VectorXd> penalty_type,
                          const Eigen::Map<Eigen::VectorXd> cmult,
                          const Eigen::Map<Eigen::VectorXd> quantiles,
                          const Rcpp::IntegerVector & num_penalty,
                          const Rcpp::NumericVector & penalty_ratio,
                          const Eigen::Map<Eigen::VectorXd> penalty_user,
                          const Eigen::Map<Eigen::VectorXd> penalty_user_ext,
                          Eigen::VectorXd lower_cl,
                          Eigen::VectorXd upper_cl,
                          const std::string & family,
                          const double & thresh,
                          const int & maxit,
                          const int & ne,
                          const int & nx) {

    // initialize objects to hold means, variances, sds of all variables
    const int n = x.rows();
    const int nv_x = x.cols();
    const int nv_fixed = fixed.size() == 0 ? 0 : fixed.cols();
    const int nv_ext = ext.size() == 0 ? 0 : ext.cols();
    const int nv_total = nv_x + nv_fixed + intr[1] + nv_ext;
    Eigen::VectorXd xm = Eigen::VectorXd::Constant(nv_total, 0.0);
    Eigen::VectorXd xv = Eigen::VectorXd::Constant(nv_total, 1.0);
    Eigen::VectorXd xs = Eigen::VectorXd::Constant(nv_total, 1.0);

    // map to correct size of fixed matrix
    Eigen::Map<const Eigen::MatrixXd> fixedmap(fixed.data(), fixed.rows(), nv_fixed);

    // scale user weights
    weights_user.array() = weights_user.array() / weights_user.sum();

    // compute moments of matrices and create XZ (if external data present)
    compute_moments(x, weights_user, xm, xv, xs, true, stnd[0], 0);
    compute_moments(fixedmap, weights_user, xm, xv, xs, true, stnd[0], nv_x);
    const Eigen::MatrixXd xz = create_XZ(x, ext, xm, xv, xs, intr[1]);
    compute_moments(xz, weights_user, xm, xv, xs, false, stnd[1], nv_x + nv_fixed);

    // choose solver based on family
    std::unique_ptr<CoordSolver<MapSpMat> > solver;
    if (family == "gaussian") {
        solver.reset(new CoordSolver<MapSpMat>(y, x, fixed, xz, xm.data(), xv.data(), xs.data(),
                                               weights_user, intr, penalty_type.data(),
                                               cmult.data(), quantiles,
                                               upper_cl.data(), lower_cl.data(),
                                               ne, nx, thresh, maxit));
    }
    else if (family == "binomial") {
        solver.reset(new BinomialSolver<MapSpMat>(y, x, fixed, xz, xm.data(), xv.data(), xs.data(),
                                                  weights_user, intr, penalty_type.data(),
                                                  cmult.data(), quantiles,
                                                  upper_cl.data(), lower_cl.data(),
                                                  ne, nx, thresh, maxit));

    }

    // Object to hold results for all penalty combinations
    const int num_combn = num_penalty[0] * num_penalty[1];
    Hierr<MapSpMat, MapMat> estimates = Hierr<MapSpMat, MapMat>(n, nv_x, nv_fixed, nv_ext, nv_total,
                                                                intr[0], intr[1], ext.data(), xm.data(),
                                                                xs.data(), num_combn);

    // compute penalty path for 1st level variables
    Eigen::VectorXd path(num_penalty[0]);

    compute_penalty(path, penalty_user, penalty_type[0],
                    penalty_ratio[0], solver->getGradient(),
                    solver->getCmult(), 0, nv_x);

    // compute penalty path for 2nd level variables
    Eigen::VectorXd path_ext(num_penalty[1]);
    if (nv_ext > 0) {
        compute_penalty(path_ext, penalty_user_ext, penalty_type[nv_x + nv_fixed + intr[1]],
                        penalty_ratio[1], solver->getGradient(),
                        solver->getCmult(), nv_x + nv_fixed + intr[1], nv_total);
    } else {
        path_ext[0] = 0.0;
    }

    // solve grid of penalties in decreasing order
    double b0_outer = solver->getBeta0();
    Eigen::VectorXd betas_outer = solver->getBetas();
    Eigen::VectorXd strong_sum = Eigen::VectorXd::Zero(num_combn);
    Eigen::VectorXd active_sum = Eigen::VectorXd::Zero(num_combn);

    int idx_pen = 0;
    for (int m = 0; m < num_penalty[0]; ++m) {
        solver->setPenalty(path[m], 0);
        for (int m2 = 0; m2 < num_penalty[1]; ++m2, ++idx_pen) {
            solver->setPenalty(path_ext[m2], 1);
            if (m2 == 0 && num_penalty[1] > 1) {
                solver->warm_start(y, b0_outer, betas_outer);
                solver->update_strong(path, path_ext, m, m2);
                solver->solve();
                b0_outer = solver->getBeta0();
                betas_outer = solver->getBetas();
            }
            else {
                solver->update_strong(path, path_ext, m, m2);
                solver->solve();
            }
            strong_sum[idx_pen] = sum(solver->getStrongSet());
            active_sum[idx_pen] = sum(solver->getActiveSet());
            estimates.add_results(solver->getBeta0(), solver->getBetas(), idx_pen);
        }
    }

    int status = 0;

    // collect results in list and return to R
    return Rcpp::List::create(Rcpp::Named("beta0") = estimates.getBeta0(),
                              Rcpp::Named("betas") = estimates.getBetas(),
                              Rcpp::Named("alpha0") = estimates.getAlpha0(),
                              Rcpp::Named("alphas") = estimates.getAlphas(),
                              Rcpp::Named("num_passes") = solver->getNumPasses(),
                              Rcpp::Named("penalty") = path,
                              Rcpp::Named("penalty_ext") = path_ext,
                              Rcpp::Named("strong_sum") = strong_sum,
                              Rcpp::Named("active_sum") = active_sum,
                              Rcpp::Named("status") = status);
}

