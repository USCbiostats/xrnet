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

#endif // HIERR_UTILS_H
