#ifndef XRNET_UTILS_H
#define XRNET_UTILS_H

#include <RcppEigen.h>
#include <string.h>

void compute_penalty(Eigen::Ref<Eigen::VectorXd> path,
                     const Eigen::Ref<const Eigen::VectorXd> & penalty_user,
                     const double & penalty_type,
                     const double & penalty_ratio,
                     const Eigen::Ref<const Eigen::VectorXd> & gradient,
                     const Eigen::Ref<const Eigen::VectorXd> & cmult,
                     const int & begin,
                     const int & end,
                     const double & ys);

double logit_inv(double x);

template <typename TX>
Eigen::MatrixXd computeResponse(const TX & X,
                                const Eigen::Ref<const Eigen::MatrixXd> & Fixed,
                                const Eigen::Ref<const Eigen::VectorXd> & beta0,
                                const Eigen::Ref<const Eigen::MatrixXd> & betas,
                                const Eigen::Ref<const Eigen::MatrixXd> & gammas,
                                const std::string & response_type,
                                const std::string & family) {

    Eigen::MatrixXd pred(X.rows(), betas.cols());
    if (gammas.cols() > 0)
        pred = Eigen::VectorXd::Constant(X.rows(), 1.0) * beta0.transpose() + X * betas + Fixed * gammas;
    else
        pred = Eigen::VectorXd::Constant(X.rows(), 1.0) * beta0.transpose() + X * betas;

    if (response_type == "response") {
        if (family == "binomial") {
            pred = pred.unaryExpr(&logit_inv);
        }
    }
    return pred;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif // XRNET_UTILS_H
