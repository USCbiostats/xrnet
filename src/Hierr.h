#ifndef HIERR_H
#define HIERR_H

#include <RcppEigen.h>
#include "CoordSolver.h"

template <typename TX, typename TZ>
class Hierr {

    typedef Eigen::VectorXd VecXd;
    typedef Eigen::MatrixXd MatXd;
    typedef Eigen::Map<const Eigen::VectorXd> MapVec;
    typedef Eigen::Map<const Eigen::MatrixXd> MapMat;

protected:
    const int n;
    const int nv_x;
    const int nv_fixed;
    const int nv_ext;
    const bool intr;
    const bool intr_ext;
    TZ ext;
    MapVec xm;
    MapVec xs;
    VecXd beta0;
    MatXd betas;
    MatXd gammas;
    VecXd alpha0;
    MatXd alphas;
    VecXd strong_sum;

public:
    // constructor (dense external)
    Hierr(const int & n_,
          const int & nv_x_,
          const int & nv_fixed_,
          const int & nv_ext_,
          const int & nv_total_,
          const bool & intr_,
          const bool & intr_ext_,
          const double * extptr,
          const double * xmptr,
          const double * xsptr,
          const int & num_penalty_) :
    n(n_),
    nv_x(nv_x_),
    nv_fixed(nv_fixed_),
    nv_ext(nv_ext_),
    intr(intr_),
    intr_ext(intr_ext_),
    ext(extptr, nv_x_, nv_ext_),
    xm(xmptr, nv_total_),
    xs(xsptr, nv_total_)
    {
        beta0 = Eigen::VectorXd::Zero(num_penalty_);
        betas = Eigen::MatrixXd::Zero(nv_x_, num_penalty_);
        gammas = Eigen::MatrixXd::Zero(nv_fixed_, num_penalty_);
        alpha0 = Eigen::VectorXd::Zero(num_penalty_);
        alphas = Eigen::MatrixXd::Zero(nv_ext_, num_penalty_);
        strong_sum = Eigen::VectorXd::Zero(num_penalty_);
    };

    // destructor
    virtual ~Hierr(){};

    // getters
    VecXd getXm(){return xm;};
    VecXd getBeta0(){return beta0;};
    double getBeta0(const int & idx){return beta0[idx];};
    MatXd getBetas(){return betas;};
    VecXd getBetas(const int & idx){return betas.col(idx);};
    VecXd getAlpha0(){return alpha0;};
    MatXd getAlphas(){return alphas;};

    // save results for single penalty
    virtual void add_results(const double & b0, VecXd coef, const int & idx) {

        // unstandardize variables
        coef = coef.cwiseProduct(xs);

        // get external coefficients
        if (nv_ext > 0) {
            alphas.col(idx) = coef.tail(nv_ext);
        }

        // unstandardize predictors w/ external data (x)
        if (nv_ext + intr_ext > 0) {
            VecXd z_alpha = Eigen::VectorXd::Zero(nv_x);
            if (intr_ext) {
                z_alpha.array() += coef[nv_x + nv_fixed];
            }
            if (nv_ext > 0) {
                z_alpha += ext * coef.tail(nv_ext);
            }
            betas.col(idx) = z_alpha * xs.head(nv_x) + coef.head(nv_x);
        }
        else {
            betas.col(idx) = coef.head(nv_x);
        }

        // unstandardize predictors w/o external data (fixed)
        if (nv_fixed > 0) {
            gammas.col(idx) = coef.segment(nv_x, nv_fixed);
        }

        // compute 2nd level intercepts
        if (intr_ext) {
            if (nv_ext > 0) {
                alpha0[idx] = coef.tail(nv_ext).sum() / nv_ext - xm.tail(nv_ext).dot(coef.tail(nv_ext));
            }
            else {
                alpha0[idx] = coef.tail(nv_ext).sum() / nv_ext;
            }
        }

        // compute 1st level intercepts
        if (intr) {
            beta0[idx] = b0 - xm.head(nv_x).dot(betas.col(idx));
            if (nv_fixed > 0) {
                beta0[idx] -= xm.segment(nv_x, nv_fixed).dot(gammas.col(idx));
            }
        }
    };
};

#endif // HIERR_H




