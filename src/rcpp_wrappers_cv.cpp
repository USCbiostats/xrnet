// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <string.h>
#include "DataFunctions.h"
#include "HierrCV.h"
#include "HierrUtils.h"
#include "CoordDescTypes.h"
#include "BinomialSolver.h"

// [[Rcpp::export]]
Eigen::VectorXd fitDenseDenseCV(SEXP x,
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
                                const std::string & user_loss,
                                const Eigen::Map<Eigen::VectorXi> test_idx,
                                const double & thresh,
                                const int & maxit,
                                const int & ne,
                                const int & nx) {

    // get pointer to big.matrix for X and map to eigen matrix
    Rcpp::XPtr<BigMatrix> xptr(x);
    Eigen::Map<const Eigen::MatrixXd> xmap((const double *)xptr->matrix(),
                                           xptr->nrow(),
                                           xptr->ncol());

    // initialize objects to hold means, variances, sds of all variables
    const int n = xmap.rows();
    const int nv_x = xmap.cols();
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
    compute_moments(xmap, weights_user, xm, xv, xs, true, stnd[0], 0);
    compute_moments(fixedmap, weights_user, xm, xv, xs, true, stnd[0], nv_x);
    const Eigen::MatrixXd xz = create_XZ(xmap, ext, xm, xv, xs, intr[1]);
    compute_moments(xz, weights_user, xm, xv, xs, false, stnd[1], nv_x + nv_fixed);

    // choose solver based on outcome
    std::unique_ptr<CoordSolver<MapMat> > solver;
    if (family == "gaussian") {
        solver.reset(new CoordSolver<MapMat>(y, xmap, fixedmap, xz, xm.data(), xv.data(), xs.data(),
                                             weights_user, intr[0], penalty_type.data(),
                                             cmult.data(), quantiles, upper_cl.data(),
                                             lower_cl.data(), ne, nx, thresh, maxit));
    }
    else if (family == "binomial") {
        solver.reset(new BinomialSolver<MapMat>(y, xmap, fixedmap, xz, xm.data(), xv.data(), xs.data(),
                                                weights_user, intr[0], penalty_type.data(),
                                                cmult.data(), quantiles, upper_cl.data(),
                                                lower_cl.data(), ne, nx, thresh, maxit));
    }

    // Object to hold results for all penalty combinations
    const int num_combn = num_penalty[0] * num_penalty[1];
    HierrCV<MapMat, MapMat> results = HierrCV<MapMat, MapMat>(n, nv_x, nv_fixed, nv_ext, nv_total,
                                                              intr[0], intr[1], ext.data(), xm.data(),
                                                              xs.data(), num_combn, family, user_loss,
                                                              test_idx, xmap, y);

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
            results.add_results(solver->getBeta0(), solver->getBetas(), idx_pen);
        }
    }

    // return results
    return results.get_error_mat();
}



// [[Rcpp::export]]
Eigen::VectorXd fitSparseDenseCV(const Eigen::MappedSparseMatrix<double> x,
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
                                 const std::string & user_loss,
                                 const Eigen::Map<Eigen::VectorXi> test_idx,
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
    HierrCV<MapSpMat, MapMat> results = HierrCV<MapSpMat, MapMat>(n, nv_x, nv_fixed, nv_ext, nv_total,
                                                                  intr[0], intr[1], ext.data(), xm.data(),
                                                                  xs.data(), num_combn, family, user_loss,
                                                                  test_idx, x, y);

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
            results.add_results(solver->getBeta0(), solver->getBetas(), idx_pen);
        }
    }

    // return error matrix
    return results.get_error_mat();
}

