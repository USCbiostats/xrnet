#ifndef GAUSSIAN_SOLVER_H
#define GAUSSIAN_SOLVER_H

#include "CoordSolver.h"

template <typename T>
class GaussianSolver : public CoordSolver<T> {

    typedef Eigen::VectorXd VecXd;
    typedef Eigen::VectorXi VecXi;
    typedef Eigen::Map<const Eigen::MatrixXd> MapMat;
    typedef Eigen::MappedSparseMatrix<double> MapSpMat;
    typedef Eigen::Map<const Eigen::VectorXd> MapVec;

private:
    using CoordSolver<T>::wgts;
    using CoordSolver<T>::wgts_user;
    using CoordSolver<T>::y;
    using CoordSolver<T>::wgts_sum;
    using CoordSolver<T>::X;
    using CoordSolver<T>::xs;
    using CoordSolver<T>::xm;
    using CoordSolver<T>::Fixed;
    using CoordSolver<T>::XZ;
    using CoordSolver<T>::gradient;
    using CoordSolver<T>::residuals;
    using CoordSolver<T>::intercept;
    using CoordSolver<T>::ym;
    using CoordSolver<T>::ys;

public:
    // constructor (dense X matrix)
    GaussianSolver(const Eigen::Ref<const Eigen::MatrixXd> & y_,
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
                   max_iterations_)
                   {
                       init();
                   };

    // constructor (sparse X matrix)
    GaussianSolver(const Eigen::Ref<const Eigen::MatrixXd> & y_,
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
                       max_iterations_)
                       {
                           init();
                       };

    // destructor
    virtual ~GaussianSolver() {}

    // initialize function
    void init() {
        wgts = wgts_user;
        wgts_sum = wgts.sum();
        ym = y.col(0).cwiseProduct(wgts_user).sum();
        ys = std::sqrt(y.col(0).cwiseProduct(y.col(0).cwiseProduct(wgts_user)).sum() - ym * ym);
        if (!intercept) {ym = 0.0;}
        residuals.array() = wgts.array() * (y.col(0).array() - ym) / ys;
        double resids_sum = residuals.sum();

        int idx = 0;
        for (int k = 0; k < X.cols(); ++k, ++idx) {
            gradient[idx] = xs[idx] * (X.col(k).dot(residuals) - xm[idx] * resids_sum);
        }
        idx += Fixed.cols();
        for (int k = 0; k < XZ.cols(); ++k, ++idx) {
            gradient[idx] = xs[idx] * (XZ.col(k).dot(residuals) - xm[idx] * resids_sum);
        }
    }
};

#endif // GAUSSIAN_SOLVER_H
