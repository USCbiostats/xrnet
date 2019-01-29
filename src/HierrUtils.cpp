#include "HierrUtils.h"
#include <bigmemory/MatrixAccessor.hpp>

typedef Eigen::Map<const Eigen::MatrixXd> MapMat;
typedef Eigen::MappedSparseMatrix<double> MapSpMat;

void compute_penalty(Eigen::Ref<Eigen::VectorXd> path,
                     const Eigen::Ref<const Eigen::VectorXd> & penalty_user,
                     const double & penalty_type,
                     const double & penalty_ratio,
                     const Eigen::Ref<const Eigen::VectorXd> & gradient,
                     const Eigen::Ref<const Eigen::VectorXd> & cmult,
                     const int & begin,
                     const int & end) {
    const int npenalty = path.size();
    if (penalty_user[0] == 0.0) {
        path[0] = 9.9e35;
        double max_penalty = 0.0;
        for (int k = begin; k < end; ++k) {
            if (cmult[k] > 0.0) {
                max_penalty = std::max(max_penalty, std::abs(gradient[k] / cmult[k]));
            }
        }
        double eqs = std::max(1e-6, penalty_ratio);
        double alf = pow(eqs, 1.0 / (npenalty - 1));
        path[1] = alf * (max_penalty / std::max(penalty_type, 0.001));
        for (int l = 2; l < npenalty; ++l) {
            path[l] = alf * path[l - 1];
        }
    }
    else {
        path = penalty_user;
    }
}

// [[Rcpp::export]]
Eigen::MatrixXd computeResponseRcpp(SEXP X,
                                    const int & mattype_x,
                                    const Eigen::Map<Eigen::MatrixXd> Fixed,
                                    const Eigen::Map<Eigen::VectorXd> beta0,
                                    const Eigen::Map<Eigen::MatrixXd> betas,
                                    const Eigen::Map<Eigen::MatrixXd> gammas) {

    if (mattype_x == 1) {
        Rcpp::NumericMatrix x_mat(X);
        MapMat xmap((const double *) &x_mat[0], x_mat.rows(), x_mat.cols());
        return computeResponse<MapMat>(xmap, Fixed, beta0, betas, gammas);
    } else if (mattype_x == 2) {
        Rcpp::XPtr<BigMatrix> xptr(X);
        MapMat xmap((const double *)xptr->matrix(), xptr->nrow(), xptr->ncol());
        return computeResponse<MapMat>(xmap, Fixed, beta0, betas, gammas);
    } else {
        return computeResponse<MapSpMat>(Rcpp::as<MapSpMat>(X), Fixed, beta0, betas, gammas);
    }
}

