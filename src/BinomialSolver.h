#ifndef BINOMIAL_SOLVER_H
#define BINOMIAL_SOLVER_H

#include "CoordSolver.h"

template <typename T>
class BinomialSolver : public CoordSolver<T> {

    typedef Eigen::VectorXd VecXd;
    typedef Eigen::VectorXi VecXi;
    typedef Eigen::Map<const Eigen::MatrixXd> MapMat;
    typedef Eigen::MappedSparseMatrix<double> MapSpMat;
    typedef Eigen::Map<const Eigen::VectorXd> MapVec;

private:
    VecXd xbeta;
    VecXd prob;
    using CoordSolver<T>::n;
    using CoordSolver<T>::nv_total;
    using CoordSolver<T>::intercept;
    using CoordSolver<T>::wgts;
    using CoordSolver<T>::wgts_user;
    using CoordSolver<T>::wgts_sum;
    using CoordSolver<T>::y;
    using CoordSolver<T>::X;
    using CoordSolver<T>::Fixed;
    using CoordSolver<T>::XZ;
    using CoordSolver<T>::xm;
    using CoordSolver<T>::xs;
    using CoordSolver<T>::xv;
    using CoordSolver<T>::residuals;
    using CoordSolver<T>::gradient;
    using CoordSolver<T>::betas;
    using CoordSolver<T>::b0;
    using CoordSolver<T>::betas_prior;
    using CoordSolver<T>::b0_prior;
    using CoordSolver<T>::tolerance_irls;
    using CoordSolver<T>::penalty;
    using CoordSolver<T>::penalty_type;
    using CoordSolver<T>::cmult;
    using CoordSolver<T>::strong_set;
    const double prob_thresh = 1e-9;
    double xbeta_thresh;


public:
    // constructor (dense X matrix)
    BinomialSolver(const Eigen::Ref<const Eigen::MatrixXd> & y_,
                   const Eigen::Ref<const Eigen::MatrixXd> & X_,
                   const Eigen::Ref<const Eigen::MatrixXd> & Fixed_,
                   const Eigen::Ref<const Eigen::MatrixXd> & XZ_,
                   const double * xmptr,
                   double * xvptr,
                   const double * xsptr,
                   VecXd wgts_user_,
                   bool intercept_,
                   const double * penalty_type_,
                   const double * cmult_,
                   VecXd quantiles_,
                   const double * ucl_,
                   const double * lcl_,
                   int ne_,
                   int nx_,
                   double tolerance_,
                   int max_iterations_) :
    CoordSolver<T>(y_,
                   X_,
                   Fixed_,
                   XZ_,
                   xmptr,
                   xvptr,
                   xsptr,
                   wgts_user_,
                   intercept_,
                   penalty_type_,
                   cmult_,
                   quantiles_,
                   ucl_,
                   lcl_,
                   ne_,
                   nx_,
                   tolerance_,
                   max_iterations_),
                   xbeta(n),
                   prob(n)
                   {
                       xbeta_thresh = log((1 - prob_thresh) / prob_thresh);
                       init();
                   };

    // constructor (sparse X matrix)
    BinomialSolver(const Eigen::Ref<const Eigen::MatrixXd> & y_,
                   const MapSpMat X_,
                   const Eigen::Ref<const Eigen::MatrixXd> & Fixed_,
                   const Eigen::Ref<const Eigen::MatrixXd> & XZ_,
                   const double * xmptr,
                   double * xvptr,
                   const double * xsptr,
                   VecXd wgts_user_,
                   bool intercept_,
                   const double * penalty_type_,
                   const double * cmult_,
                   VecXd quantiles_,
                   const double * ucl_,
                   const double * lcl_,
                   int ne_,
                   int nx_,
                   double tolerance_,
                   int max_iterations_) :
        CoordSolver<T>(y_,
                       X_,
                       Fixed_,
                       XZ_,
                       xmptr,
                       xvptr,
                       xsptr,
                       wgts_user_,
                       intercept_,
                       penalty_type_,
                       cmult_,
                       quantiles_,
                       ucl_,
                       lcl_,
                       ne_,
                       nx_,
                       tolerance_,
                       max_iterations_),
                       xbeta(n),
                       prob(n)
                       {
                           xbeta_thresh = log((1 - prob_thresh) / prob_thresh);
                           init();
                       };

    // destructor
    virtual ~BinomialSolver() {}

    // initialize function
    void init() {

        xbeta = Eigen::VectorXd::Zero(n);
        prob = Eigen::VectorXd::Constant(n, 0.5);

        // initial p(y = 1) for all obs.
        double prob0;
        if (intercept) {
            prob0 = wgts_user.dot(y.col(0));
        } else {
            prob0 = 0.5;
        }

        // initial value of intercept and sum of wgts
        b0 = log(prob0 / (1.0 - prob0));
        wgts_sum = prob0 * (1 - prob0);

        // initial wgts
        wgts.array() = wgts_user.array() * wgts_sum;

        // initial residuals
        residuals.array() = wgts_user.array() * (y.col(0).array() - prob0);

        // initial weighted sum squares x / xz cols and gradient
        int idx = 0;
        for (int k = 0; k < X.cols(); ++k, ++idx) {
            gradient[idx] = xs[idx] * (X.col(k).dot(residuals) - xm[idx] * residuals.sum());
            xv[idx] = std::pow(xs[idx], 2) * (X.col(k).cwiseProduct(X.col(k)) - 2 * xm[idx] * X.col(k) + std::pow(xm[idx], 2) * Eigen::VectorXd::Ones(n)).adjoint() * wgts;
        }
        for (int k = 0; k < Fixed.cols(); ++k, ++idx) {
            gradient[idx] = xs[idx] * (Fixed.col(k).dot(residuals) - xm[idx] * residuals.sum());
            xv[idx] = std::pow(xs[idx], 2) * (Fixed.col(k).cwiseProduct(Fixed.col(k)) - 2 * xm[idx] * Fixed.col(k) + std::pow(xm[idx], 2) * Eigen::VectorXd::Ones(n)).adjoint() * wgts;
        }
        for (int k = 0; k < XZ.cols(); ++k, ++idx) {
            gradient[idx] = xs[idx] * (XZ.col(k).dot(residuals) - xm[idx] * residuals.sum());
            xv[idx] = std::pow(xs[idx], 2) * (XZ.col(k).cwiseProduct(XZ.col(k)) - 2 * xm[idx] * XZ.col(k) + std::pow(xm[idx], 2) * Eigen::VectorXd::Ones(n)).adjoint() * wgts;
        }
    }

    // warm start initialization given current estimates
    virtual void warm_start(const double & b0_start,
                            const Eigen::Ref<const Eigen::VectorXd> & betas_start) {

        // initialize estimates with provided values
        b0 = b0_start;
        betas = betas_start;

        // update residuals, working response, weighted sum squres X / XZ
        update_quadratic();

        // update gradients given current residuals
        int idx = 0;
        for (int k = 0; k < X.cols(); ++k, ++idx) {
            gradient[idx] = xs[idx] * (X.col(k).dot(residuals) - xm[idx] * residuals.sum());
        }
        idx += Fixed.cols();
        for (int k = 0; k < XZ.cols(); ++k, ++idx) {
            gradient[idx] = xs[idx] * (XZ.col(k).dot(residuals) - xm[idx] * residuals.sum());
        }
    }

    // update quadratic approx. of log-likelihood
    virtual void update_quadratic() {
        // compute linear predictor (X * beta)
        xbeta.array() = b0;
        int idx = 0;
        for (int j = 0; j < X.cols(); ++j, ++idx) {
            if (strong_set[idx] && betas[idx] != 0.0) {
                xbeta += xs[idx] * (X.col(j) - xm[idx] * Eigen::VectorXd::Ones(n)) * betas[idx];
            }
        }
        for (int j = 0; j < Fixed.cols(); ++j, ++idx) {
            xbeta += xs[idx] * (Fixed.col(j) - xm[idx] * Eigen::VectorXd::Ones(n)) * betas[idx];
        }
        for (int j = 0; j < XZ.cols(); ++j, ++idx) {
            if (strong_set[idx] && betas[idx] != 0.0) {
                xbeta += xs[idx] * (XZ.col(j) - xm[idx] * Eigen::VectorXd::Ones(n)) * betas[idx];
            }
        }

        // compute predicted probabilities
        for (int i = 0; i < n; ++i) {
            prob[i] = abs(xbeta[i]) < xbeta_thresh ? 1.0 / (1.0 + exp(-xbeta[i])) : xbeta[i] > 0 ? 1.0 : 0.0;
        }

        // update weights
        wgts.array() = wgts_user.array() * prob.array() * (1 - prob.array());
        wgts_sum = wgts.sum();

        // update residuals
        residuals.array() = wgts_user.array() * (y.col(0).array() - prob.array());

        // update weighted sum squares x / xz cols
        idx = 0;
        for (int k = 0; k < X.cols(); ++k, ++idx) {
            if (strong_set[idx])
                xv[idx] = std::pow(xs[idx], 2) * (X.col(k).cwiseProduct(X.col(k)) - 2 * xm[idx] * X.col(k) + std::pow(xm[idx], 2) * Eigen::VectorXd::Ones(n)).adjoint() * wgts;
        }
        for (int k = 0; k < Fixed.cols(); ++k, ++idx) {
            if (strong_set[idx])
                xv[idx] = std::pow(xs[idx], 2) * (Fixed.col(k).cwiseProduct(Fixed.col(k)) - 2 * xm[idx] * Fixed.col(k) + std::pow(xm[idx], 2) * Eigen::VectorXd::Ones(n)).adjoint() * wgts;
        }
        for (int k = 0; k < XZ.cols(); ++k, ++idx) {
            if (strong_set[idx])
                xv[idx] = std::pow(xs[idx], 2) * (XZ.col(k).cwiseProduct(XZ.col(k)) - 2 * xm[idx] * XZ.col(k) + std::pow(xm[idx], 2) * Eigen::VectorXd::Ones(n)).adjoint() * wgts;
        }
    }

    // check convergence of IRLS
    virtual bool converged() {
        bool converged_outer = true;
        if (wgts.sum() >= prob_thresh) {
            if (wgts.sum() * std::pow(b0 - b0_prior, 2) > tolerance_irls) {
                converged_outer = false;
            }
            else {
                for (int k = 0; k < nv_total; ++k) {
                    if (strong_set[k] && xv[k] * std::pow(betas[k] - betas_prior[k], 2) > tolerance_irls) {
                        converged_outer = false;
                        break;
                    }
                }
            }
        }
        betas_prior = betas;
        b0_prior = b0;
        return converged_outer;
    }

    // check kkt conditions
    virtual bool check_kkt() {
        int num_violations = 0;
        int idx = 0;
        double resid_sum = residuals.sum();
        for (int k = 0; k < X.cols(); ++k, ++idx) {
            if (!strong_set[idx]) {
                gradient[idx] = xs[idx] * (X.col(k).dot(residuals) - xm[idx] * resid_sum);
                if (std::abs(gradient[idx]) > penalty[0] * penalty_type[idx] * cmult[idx]) {
                    strong_set[idx] = true;
                    xv[idx] = std::pow(xs[idx], 2) * (X.col(k).cwiseProduct(X.col(k)) - 2 * xm[idx] * X.col(k) + std::pow(xm[idx], 2) * Eigen::VectorXd::Ones(n)).adjoint() * wgts;
                    ++num_violations;
                }
            }
        }
        idx += Fixed.cols();
        for (int k = 0; k < XZ.cols(); ++k, ++idx) {
            if (!strong_set[idx]) {
                gradient[idx] = xs[idx] * (XZ.col(k).dot(residuals) - xm[idx] * resid_sum);
                if (std::abs(gradient[idx]) > penalty[1] * penalty_type[idx] * cmult[idx]) {
                    strong_set[idx] = true;
                    xv[idx] = std::pow(xs[idx], 2) * (XZ.col(k).cwiseProduct(XZ.col(k)) - 2 * xm[idx] * XZ.col(k) + std::pow(xm[idx], 2) * Eigen::VectorXd::Ones(n)).adjoint() * wgts;
                    ++num_violations;
                }
            }
        }
        return num_violations == 0;
    }
};

#endif // BINOMIAL_SOLVER_H
