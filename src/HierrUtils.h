#ifndef HIERR_UTILS_H
#define HIERR_UTILS_H

#include <RcppEigen.h>

void compute_penalty(Eigen::Ref<Eigen::VectorXd> path,
                     const Eigen::Ref<const Eigen::VectorXd> & penalty_user,
                     const double & penalty_type,
                     const double & penalty_ratio,
                     const Eigen::Ref<const Eigen::VectorXd> & gradient,
                     const Eigen::Ref<const Eigen::VectorXd> & cmult,
                     const int & begin,
                     const int & end);

template <typename TX>
Eigen::MatrixXd computeResponse(const TX & X,
                                const Eigen::Ref<const Eigen::MatrixXd> & Fixed,
                                const Eigen::Ref<const Eigen::VectorXd> & beta0,
                                const Eigen::Ref<const Eigen::MatrixXd> & betas,
                                const Eigen::Ref<const Eigen::MatrixXd> & gammas) {

    if (gammas.cols() > 0)
        return Eigen::VectorXd::Constant(X.rows(), 1.0) * beta0.transpose() + X * betas + Fixed * gammas;
    else
        return Eigen::VectorXd::Constant(X.rows(), 1.0) * beta0.transpose() + X * betas;
}

#endif // HIERR_UTILS_H
