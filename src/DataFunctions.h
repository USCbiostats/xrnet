#ifndef DATA_FUNCTIONS_H
#define DATA_FUNCTIONS_H

#include <RcppEigen.h>

template <typename matType>
void compute_moments(const matType & X,
                     const Eigen::Ref<const Eigen::VectorXd> & wgts_user,
                     Eigen::Ref<Eigen::VectorXd> xm,
                     Eigen::Ref<Eigen::VectorXd> cent,
                     Eigen::Ref<Eigen::VectorXd> xv,
                     Eigen::Ref<Eigen::VectorXd> xs,
                     const bool & centered,
                     const bool & scaled,
                     int idx) {
    if (centered) {
        if (scaled) {
            for (int j = 0; j < X.cols(); ++j, ++idx) {
                auto xj = X.col(j);
                xm[idx] = xj.cwiseProduct(wgts_user).sum();
                cent[idx] = xm[idx];
                xs[idx] = 1 / std::sqrt(xj.cwiseProduct(xj.cwiseProduct(wgts_user)).sum() - xm[idx] * xm[idx]);
            }
        } else {
            for (int j = 0; j < X.cols(); ++j, ++idx) {
                auto xj = X.col(j);
                xm[idx] = xj.cwiseProduct(wgts_user).sum();
                cent[idx] = xm[idx];
                xv[idx] = xj.cwiseProduct(xj.cwiseProduct(wgts_user)).sum() - xm[idx] * xm[idx];
            }
        }
    }
    else {
        if (scaled) {
            for (int j = 0; j < X.cols(); ++j, ++idx) {
                auto xj = X.col(j);
                xm[idx] = xj.cwiseProduct(wgts_user).sum();
                double vc = xj.cwiseProduct(xj.cwiseProduct(wgts_user)).sum() - xm[idx] * xm[idx];
                xs[idx] = 1 / std::sqrt(vc);
                xv[idx] = 1.0 + xm[idx] * xm[idx] / vc;
            }
        } else {
            for (int j = 0; j < X.cols(); ++j, ++idx) {
                auto xj = X.col(j);
                xm[idx] = xj.cwiseProduct(wgts_user).sum();
                xv[idx] = xj.cwiseProduct(xj.cwiseProduct(wgts_user)).sum();
            }
        }
    }
}

template <typename matA, typename matB>
Eigen::MatrixXd create_XZ(const matA & X,
                          const matB & Z,
                          Eigen::Ref<Eigen::VectorXd> xm,
                          const Eigen::Ref<const Eigen::VectorXd> & cent,
                          const Eigen::Ref<const Eigen::VectorXd> & wgts_user,
                          Eigen::Ref<Eigen::VectorXd> xv,
                          Eigen::Ref<Eigen::VectorXd> xs,
                          const bool & intr_ext,
                          const bool & scale_z,
                          int idx) {

    // initialize XZ matrix
    Eigen::MatrixXd XZ(0, 0);

    if (Z.size() == 0)
        return XZ;
    else
        XZ.resize(X.rows(), intr_ext + Z.cols());

    // map means and sds of X
    auto cent_x = cent.head(X.cols());
    auto xs_x = xs.head(X.cols());

    int col_xz = 0;

    // add intercept
    if (intr_ext) {
        auto xzj = XZ.col(col_xz);
        xzj = (X * xs_x).array() - xs_x.cwiseProduct(cent_x).sum();
        double xzj_mean = xzj.cwiseProduct(wgts_user).sum();
        xv[idx] = xzj.cwiseProduct(xzj.cwiseProduct(wgts_user)).sum() - xzj_mean * xzj_mean;
        ++idx;
        ++col_xz;
    }

    // fill in columns of XZ
    for (int j = 0; j < Z.cols(); ++j, ++col_xz, ++idx) {
        auto zj = Z.col(j);
        auto xzj = XZ.col(col_xz);
        xm[idx] = zj.sum() / zj.size();
        if (scale_z) {
            xs[idx] = 1 / std::sqrt(zj.cwiseProduct(zj / zj.size()).sum() - xm[idx] * xm[idx]);
        }
        xzj = xs[idx] * ((X * zj.cwiseProduct(xs_x)) - xs_x.cwiseProduct(cent_x.cwiseProduct(zj)).sum() * Eigen::VectorXd::Ones(X.rows()));
        double xzj_mean = xzj.cwiseProduct(wgts_user).sum();
        xv[idx] = xzj.cwiseProduct(xzj.cwiseProduct(wgts_user)).sum() - xzj_mean * xzj_mean;
        xzj /= xs[idx];

    }
    return XZ;
}

#endif // DATA_FUNCTIONS_H
