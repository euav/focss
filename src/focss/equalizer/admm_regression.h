#ifndef FOCSS_EQUALIZER_ADMM_REGRESSION_H_
#define FOCSS_EQUALIZER_ADMM_REGRESSION_H_

#include <armadillo>

namespace focss {
arma::cx_vec cx_soft_threshold(const arma::cx_vec& values,
                               const double& threshold);

arma::cx_vec cx_lasso_admm(const arma::cx_mat& A,
                           const arma::cx_vec& b,
                           const double& lambda);

arma::cx_vec cx_lad_admm(const arma::cx_mat& A, const arma::cx_vec& b);
}  // namespace focss

#endif  // FOCSS_EQUALIZER_ADMM_REGRESSION_H_
