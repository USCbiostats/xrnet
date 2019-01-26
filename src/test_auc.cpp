#include <RcppEigen.h>

// [[Rcpp::export]]
double auc_test(const Eigen::Map<Eigen::VectorXd> & actual,
                const Eigen::Map<Eigen::VectorXd> & pred,
                const Eigen::Map<Eigen::VectorXi> & idx) {

    int test_size = idx.size();
    Eigen::VectorXd actual_sub(test_size);
    Eigen::VectorXd pred_sub(test_size);
    for (int i = 0; i < test_size; ++i) {
        actual_sub[i] = actual[idx[i]];
        pred_sub[i] = pred[idx[i]];
    }

    const int n = pred_sub.size();
    std::vector<size_t> indx(n);
    std::iota(indx.begin(), indx.end(), 0);
    std::sort(indx.begin(), indx.end(), [&pred_sub](int i1, int i2) {return pred_sub[i1] < pred_sub[i2];});

    int n1 = 0;
    double rank_sum = 0;
    for (size_t i = 0; i < n; ++i) {
        if (actual_sub[indx[i]] == 1) {
            ++n1;
            rank_sum += i + 1;
        }
    }
    double u_value = rank_sum - (n1 * (n1 + 1)) / 2.0;
    return u_value / (n1 * (n - n1));
}

/*** R
data(BinomialExample)
fit <- glmnet(x, y, family = "binomial")
pred <- drop(predict(fit, s = fit$lambda[10], newx = x, type = "response"))
test_idx <- sample(1:length(y), 20, replace = F)
test_idx2 <- as.integer(test_idx - 1)

pROC::auc(y[test_idx], pred[test_idx])
hierr:::auc_test(as.numeric(y), pred, test_idx2)


microbenchmark::microbenchmark(
    pROC::auc(y, pred),
    hierr:::auc_test(as.numeric(y), pred, test_idx2)
)

*/
