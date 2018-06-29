#ifndef FOCSS_EQUALIZER_ADMM_REGRESSION_H_
#define FOCSS_EQUALIZER_ADMM_REGRESSION_H_

#include <armadillo>

namespace focss {
arma::cx_vec soft_threshold(const arma::cx_vec& values,
                            const double& threshold);

void soft_inplace(const double& threshold, arma::cx_vec& values);

arma::cx_vec cx_lasso_admm(const arma::cx_mat& A,
                           const arma::cx_vec& b,
                           const double& lambda);

arma::cx_vec cx_lad_admm(const arma::cx_mat& A, const arma::cx_vec& b);

arma::cx_vec cx_ridge_cholesky(const arma::cx_mat& A,
                               const arma::cx_vec& b,
                               const double& alpha);

arma::cx_vec cx_modulus_cholesky(const arma::cx_mat& A,
                                 const arma::cx_vec& b,
                                 const double& needed);

arma::cx_vec cx_ridge_qr(const arma::cx_mat& A,
                         const arma::cx_vec& b,
                         const double& alpha);

arma::cx_vec cx_sag(const arma::cx_mat& A,
                    const arma::cx_vec& b,
                    const double& lambda);
}  // namespace focss

#endif  // FOCSS_EQUALIZER_ADMM_REGRESSION_H_
