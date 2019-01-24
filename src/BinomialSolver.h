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

public:
    // constructor (dense X matrix)
    BinomialSolver(const Eigen::Ref<const Eigen::VectorXd> & y_,
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
                   xbeta(this->n),
                   prob(this->n)
                   {
                       xbeta = Eigen::VectorXd::Zero(this->n);
                       prob = Eigen::VectorXd::Constant(this->n, 0.5);
                       init();
                   };

    // constructor (sparse X matrix)
    BinomialSolver(const Eigen::Ref<const Eigen::VectorXd> & y_,
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
                       xbeta(this->n),
                       prob(this->n)
                       {
                           xbeta = Eigen::VectorXd::Zero(this->n);
                           prob = Eigen::VectorXd::Constant(this->n, 0.5);
                           init();
                       };

    // destructor
    virtual ~BinomialSolver() {};

    // initialize function
    virtual void init() {
        // initial p(y = 1) for all obs.
        double prob0;
        if (this->intercept) {
            prob0 = this->wgts_user.dot(this->y);
        } else {
            prob0 = 0.5;
        }

        // initial value of intercept and sum of wgts
        this->b0 = log(prob0 / (1.0 - prob0));
        this->wgts_sum = prob0 * (1 - prob0);

        // initial wgts
        this->wgts.array() = this->wgts_user.array() * this->wgts_sum;

        // initial residuals
        this->residuals.array() = this->wgts_user.array() * (this->y.array() - prob0);

        // initial weighted sum squares x / xz cols and gradient
        int idx = 0;
        for (int k = 0; k < this->X.cols(); ++k, ++idx) {
            this->gradient[idx] = this->xs[idx] * (this->X.col(k).dot(this->residuals) - this->xm[idx] * this->residuals.sum());
            this->xv[idx] = std::pow(this->xs[idx], 2) * (this->X.col(k).cwiseProduct(this->X.col(k)) - 2 * this->xm[idx] * this->X.col(k) + std::pow(this->xm[idx], 2) * Eigen::VectorXd::Ones(this->n)).adjoint() * this->wgts;
        }
        for (int k = 0; k < this->Fixed.cols(); ++k, ++idx) {
            this->gradient[idx] = this->xs[idx] * (this->Fixed.col(k).dot(this->residuals) - this->xm[idx] * this->residuals.sum());
            this->xv[idx] = std::pow(this->xs[idx], 2) * (this->Fixed.col(k).cwiseProduct(this->Fixed.col(k)) - 2 * this->xm[idx] * this->Fixed.col(k) + std::pow(this->xm[idx], 2) * Eigen::VectorXd::Ones(this->n)).adjoint() * this->wgts;
        }
        for (int k = 0; k < this->XZ.cols(); ++k, ++idx) {
            this->gradient[idx] = this->xs[idx] * (this->XZ.col(k).dot(this->residuals) - this->xm[idx] * this->residuals.sum());
            this->xv[idx] = std::pow(this->xs[idx], 2) * (this->XZ.col(k).cwiseProduct(this->XZ.col(k)) - 2 * this->xm[idx] * this->XZ.col(k) + std::pow(this->xm[idx], 2) * Eigen::VectorXd::Ones(this->n)).adjoint() * this->wgts;
        }
    };

    // warm start initialization given current estimates
    virtual void warm_start(const Eigen::Ref<const Eigen::VectorXd> & y,
                            const double & b0_start,
                            const Eigen::Ref<const Eigen::VectorXd> & betas_start) {

        // initialize estimates with provided values
        this->b0 = b0_start;
        this->betas = betas_start;

        // update residuals, working response, weighted sum squres X / XZ
        update_quadratic();

        // update partial gradients given current residuals
        int idx = 0;
        for (int k = 0; k < this->X.cols(); ++k, ++idx) {
            this->gradient[idx] = this->xs[idx] * (this->X.col(k).dot(this->residuals) - this->xm[idx] * this->residuals.sum());
        }
        idx += this->Fixed.cols();
        for (int k = 0; k < this->XZ.cols(); ++k, ++idx) {
            this->gradient[idx] = this->xs[idx] * (this->XZ.col(k).dot(this->residuals) - this->xm[idx] * this->residuals.sum());
        }
    };

    // update quadratic approx. of log-likelihood
    virtual void update_quadratic() {
        // compute linear predictor (X * beta)
        xbeta.array() = this->b0;
        int idx = 0;
        for (int j = 0; j < this->X.cols(); ++j, ++idx) {
            if (this->strong_set[idx]) {
                xbeta += this->xs[idx] * (this->X.col(j) - this->xm[idx] * Eigen::VectorXd::Ones(this->n)) * this->betas[idx];
            }
        }
        for (int j = 0; j < this->Fixed.cols(); ++j, ++idx) {
            xbeta += this->xs[idx] * (this->Fixed.col(j) - this->xm[idx] * Eigen::VectorXd::Ones(this->n)) * this->betas[idx];
        }
        for (int j = 0; j < this->XZ.cols(); ++j, ++idx) {
            if (this->strong_set[idx]) {
                xbeta += this->xs[idx] * (this->XZ.col(j) - this->xm[idx] * Eigen::VectorXd::Ones(this->n)) * this->betas[idx];
            }
        }

        // compute predicted probabilities
        for (int i = 0; i < this->n; ++i) {
            if (xbeta[i] > 20.72327) {
                prob[i] = 1.0;
            } else if (xbeta[i] < -20.72327) {
                prob[i] = 0.0;
            } else {
                prob[i] = 1.0 / (1.0 + exp(-xbeta[i])) ;
            }
        }

        // update weights
        this->wgts.array() = this->wgts_user.array() * prob.array() * (1 - prob.array());
        this->wgts_sum = this->wgts.sum();

        // update residuals
        this->residuals.array() = this->wgts_user.array() * (this->y.array() - prob.array());

        // update weighted sum squares x / xz cols
        idx = 0;
        for (int k = 0; k < this->X.cols(); ++k, ++idx) {
            if (this->strong_set[idx])
                this->xv[idx] = std::pow(this->xs[idx], 2) * (this->X.col(k).cwiseProduct(this->X.col(k)) - 2 * this->xm[idx] * this->X.col(k) + std::pow(this->xm[idx], 2) * Eigen::VectorXd::Ones(this->n)).adjoint() * this->wgts;
        }
        for (int k = 0; k < this->Fixed.cols(); ++k, ++idx) {
            if (this->strong_set[idx])
                this->xv[idx] = std::pow(this->xs[idx], 2) * (this->Fixed.col(k).cwiseProduct(this->Fixed.col(k)) - 2 * this->xm[idx] * this->Fixed.col(k) + std::pow(this->xm[idx], 2) * Eigen::VectorXd::Ones(this->n)).adjoint() * this->wgts;
        }
        for (int k = 0; k < this->XZ.cols(); ++k, ++idx) {
            if (this->strong_set[idx])
                this->xv[idx] = std::pow(this->xs[idx], 2) * (this->XZ.col(k).cwiseProduct(this->XZ.col(k)) - 2 * this->xm[idx] * this->XZ.col(k) + std::pow(this->xm[idx], 2) * Eigen::VectorXd::Ones(this->n)).adjoint() * this->wgts;
        }
    };

    // check convergence of IRLS
    virtual bool converged() {
        bool converged_outer = true;
        if (this->wgts.sum() >= 1e-9) {
            if (this->wgts.sum() * std::pow(this->b0 - this->b0_prior, 2) > this->tolerance_irls) {
                converged_outer = false;
            }
            else {
                for (int k = 0; k < this->nv_total; ++k) {
                    if (this->strong_set[k] && this->xv[k] * std::pow(this->betas[k] - this->betas_prior[k], 2) > this->tolerance_irls) {
                        converged_outer = false;
                        break;
                    }
                }
            }
        }
        this->betas_prior = this->betas;
        this->b0_prior = this->b0;
        return converged_outer;
    };

    // check kkt conditions
    virtual bool check_kkt() {
        int num_violations = 0;
        int idx = 0;
        double resid_sum = this->residuals.sum();
        for (int k = 0; k < this->X.cols(); ++k, ++idx) {
            if (!this->strong_set[idx]) {
                this->gradient[idx] = this->xs[idx] * (this->X.col(k).dot(this->residuals) - this->xm[idx] * resid_sum);
                if (std::abs(this->gradient[idx]) > this->penalty[0] * this->penalty_type[idx] * this->cmult[idx]) {
                    this->strong_set[idx] = true;
                    this->xv[idx] = std::pow(this->xs[idx], 2) * (this->X.col(k).cwiseProduct(this->X.col(k)) - 2 * this->xm[idx] * this->X.col(k) + std::pow(this->xm[idx], 2) * Eigen::VectorXd::Ones(this->n)).adjoint() * this->wgts;
                    ++num_violations;
                }
            }
        }
        idx += this->Fixed.cols();
        for (int k = 0; k < this->XZ.cols(); ++k, ++idx) {
            if (!this->strong_set[idx]) {
                this->gradient[idx] = this->xs[idx] * (this->XZ.col(k).dot(this->residuals) - this->xm[idx] * resid_sum);
                if (std::abs(this->gradient[idx]) > this->penalty[0] * this->penalty_type[idx] * this->cmult[idx]) {
                    this->strong_set[idx] = true;
                    this->xv[idx] = std::pow(this->xs[idx], 2) * (this->XZ.col(k).cwiseProduct(this->XZ.col(k)) - 2 * this->xm[idx] * this->XZ.col(k) + std::pow(this->xm[idx], 2) * Eigen::VectorXd::Ones(this->n)).adjoint() * this->wgts;
                    ++num_violations;
                }
            }
        }
        return num_violations == 0;
    };
};

#endif // BINOMIAL_SOLVER_H
