#ifndef ADMM_REGRESSION_H_
#define ADMM_REGRESSION_H_

#include <armadillo>

arma::cx_vec cx_soft_threshold(const arma::cx_vec& values,
                               const double& threshold);

arma::cx_vec cx_lasso_admm(const arma::cx_mat& A,
                           const arma::cx_vec& b,
                           const double& lambda);

arma::cx_vec cx_lad_admm(const arma::cx_mat& A, const arma::cx_vec& b);

#endif  // ADMM_REGRESSION_H_
