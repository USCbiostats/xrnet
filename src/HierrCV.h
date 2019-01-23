#ifndef HIERR_CV_H
#define HIERR_CV_H

#include <RcppEigen.h>
#include <unordered_map>
#include "Hierr.h"

template <typename TX, typename TZ>
class HierrCV : public Hierr<TX, TZ>  {

    typedef Eigen::VectorXd VecXd;
    typedef Eigen::MatrixXd MatXd;
    typedef Eigen::Map<const Eigen::VectorXd> MapVec;
    typedef Eigen::Map<const Eigen::MatrixXd> MapMat;
    typedef Eigen::MappedSparseMatrix<double> MapSpMat;
    typedef double (*lossPtr)(const Eigen::Ref<const Eigen::VectorXd> &,
                              const Eigen::Ref<const Eigen::VectorXd> &,
                              const Eigen::Ref<const Eigen::VectorXi> &);

protected:
    Eigen::Map<const Eigen::VectorXi> test_idx;
    TX X;
    MapVec y;
    VecXd error_mat;
    lossPtr loss_func;

public:
    // constructor (dense X)
    HierrCV(const int & n_,
            const int & nv_x_,
            const int & nv_fixed_,
            const int & nv_ext_,
            const int & nv_total_,
            const bool & intr_,
            const bool & intr_ext_,
            const TZ ext_,
            const double * xmptr,
            const double * xsptr,
            const int & num_penalty_,
            const std::string & family_,
            const std::string & user_loss_,
            const Eigen::Ref<const Eigen::VectorXi> & test_idx_,
            const Eigen::Ref<const Eigen::MatrixXd> & X_,
            const Eigen::Ref<const Eigen::VectorXd> & y_) :
    Hierr<TX, TZ>(
            n_,
            nv_x_,
            nv_fixed_,
            nv_ext_,
            nv_total_,
            intr_,
            intr_ext_,
            ext_,
            xmptr,
            xsptr,
            1),
            test_idx(test_idx_.data(), test_idx_.size()),
            X(X_.data(), n_, X_.cols()),
            y(y_.data(), n_)
            {
                error_mat = Eigen::VectorXd::Zero(num_penalty_);
                loss_func = select_loss(family_, user_loss_);
            };

    // constructor (sparse X)
    HierrCV(const int & n_,
            const int & nv_x_,
            const int & nv_fixed_,
            const int & nv_ext_,
            const int & nv_total_,
            const bool & intr_,
            const bool & intr_ext_,
            const TZ ext_,
            const double * xmptr,
            const double * xsptr,
            const int & num_penalty_,
            const std::string & family_,
            const std::string & user_loss_,
            const Eigen::Ref<const Eigen::VectorXi> & test_idx_,
            const MapSpMat X_,
            const Eigen::Ref<const Eigen::VectorXd> & y_) :
        Hierr<TX, TZ>(
                n_,
                nv_x_,
                nv_fixed_,
                nv_ext_,
                nv_total_,
                intr_,
                intr_ext_,
                ext_,
                xmptr,
                xsptr,
                1),
                test_idx(test_idx_.data(), test_idx_.size()),
                X(X_),
                y(y_.data(), n_)
                {
                    error_mat = Eigen::VectorXd::Zero(num_penalty_);
                    loss_func = select_loss(family_, user_loss_);
                };

    // destructor
    virtual ~HierrCV(){};

    // getters
    MatXd get_error_mat(){return error_mat;};

    // save results for single penalty
    virtual void add_results(const double & b0, VecXd coef, const int & idx) {

        // unstandardize variables
        coef = coef.cwiseProduct(this->xs);

        // get external coefficients
        if (this->nv_ext > 0) {
            this->alphas.col(0) = coef.tail(this->nv_ext);
        }

        // unstandardize predictors w/ external data (x)
        if (this->nv_ext + this->intr_ext > 0) {
            VecXd z_alpha = Eigen::VectorXd::Zero(this->nv_x);
            if (this->intr_ext) {
                z_alpha.array() += coef[this->nv_x + this->nv_fixed];
            }
            if (this->nv_ext > 0) {
                z_alpha += this->ext * coef.tail(this->nv_ext);
            }
            this->betas.col(0) = z_alpha.cwiseProduct(this->xs.head(this->nv_x)) + coef.head(this->nv_x);
        }
        else {
            this->betas.col(0) = coef.head(this->nv_x);
        }

        // unstandardize predictors w/o external data (fixed)
        if (this->nv_fixed > 0) {
            this->gammas.col(0) = coef.segment(this->nv_x, this->nv_fixed);
        }

        // compute 1st level intercept
        if (this->intr) {
            this->beta0[0] = b0 - this->xm.head(this->nv_x).dot(this->betas.col(0));
            if (this->nv_fixed > 0) {
                this->beta0[0] -= this->xm.segment(this->nv_x, this->nv_fixed).dot(this->gammas.col(0));
            }
        }

        // compute predicted values
        VecXd yhat = this->beta0[0] + (X * this->betas.col(0)).array();

        // compute error for test data
        error_mat[idx] = loss_func(y, yhat, test_idx);
    };

    // returns pointer to appropriate loss function based on family / user_loss
    lossPtr select_loss(const std::string & family, const std::string & user_loss) {

        /* map of string name - function name
        std::unordered_map<std::string, lossPtr> lossMap = {
            {"mse", mean_squared_error},
            {"mae", mean_absolute_error}
        };
        */

        std::unordered_map<std::string, lossPtr> lossMap;
        lossMap.insert(std::make_pair<std::string, lossPtr>("mse", mean_squared_error));
        lossMap.insert(std::make_pair<std::string, lossPtr>("mae", mean_absolute_error));

        lossPtr loss_func = nullptr;
        if (user_loss == "default")
        {
            if (family == "gaussian")
                loss_func =  mean_squared_error;
            else if (family == "binomial")
                loss_func = mean_squared_error;
        }
        else
        {
            std::vector<std::string> loss_options;
            if (family == "gaussian") {
                Rcpp::StringVector loss_options;
                loss_options.push_back("mse");
                loss_options.push_back("mae");
            }
            auto it = std::find(loss_options.begin(), loss_options.end(), user_loss);
            loss_func = lossMap[*it];
        }
        return loss_func;
    }

    static double mean_squared_error(const Eigen::Ref<const Eigen::VectorXd> & actual,
                                     const Eigen::Ref<const Eigen::VectorXd> & predicted,
                                     const Eigen::Ref<const Eigen::VectorXi> & test_idx) {
        double error = 0.0;
        for (int i = 0; i < test_idx.size(); ++i) {
            error += std::pow(actual[i] - predicted[i], 2);
        }
        return error / test_idx.size();
    }

    static double mean_absolute_error(const Eigen::Ref<const Eigen::VectorXd> & actual,
                                      const Eigen::Ref<const Eigen::VectorXd> & predicted,
                                      const Eigen::Ref<const Eigen::VectorXi> & test_idx) {
        double error = 0.0;
        for (int i = 0; i < test_idx.size(); ++i) {
            error += std::abs(actual[i] - predicted[i]);
        }
        return error / test_idx.size();
    }
};

#endif // HIERR_CV_H




