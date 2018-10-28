#ifndef FOCSS_EQUALIZER_WDM_SINO_H_
#define FOCSS_EQUALIZER_WDM_SINO_H_

#include <armadillo>
#include "focss/field.h"

namespace focss {
class WdmPPE {
  public:
    WdmPPE();
    WdmPPE(const int& radii);
    WdmPPE(const int& c_radius, const int& d_radius);
    WdmPPE& set_c_radius(const int& c_radius);
    WdmPPE& set_d_radius(const int& d_radius);

    void train(const arma::cx_mat& tx, const arma::cx_mat& rx);

    complex_t get_intercept() const;
    complex_t get_linear_weight() const;
    complex_t get_coo_weight() const;
    arma::cx_mat linear(const arma::cx_mat& field) const;
    arma::cx_mat spm(const arma::cx_mat& field) const;
    arma::cx_mat xpm(const arma::cx_mat& field) const;
    arma::cx_mat cubic(const arma::cx_mat& field) const;
    arma::cx_mat model(const arma::cx_mat& field) const;
    arma::cx_mat equalize(const arma::cx_mat& tx, const arma::cx_mat& rx) const;
    arma::cx_mat equalize(const arma::cx_mat& tx_spm,
                          const arma::cx_mat& tx_xpm,
                          const arma::cx_mat& rx) const;

  private:
    static complex_t TT(const arma::cx_mat& field,
                        const int& pol,
                        const int& k,
                        const int& m,
                        const int& n);

    static complex_t T(const arma::cx_mat& field,
                       const int& pol,
                       const int& k,
                       const int& m,
                       const int& n);

    static complex_t T1(const arma::cx_mat& field,
                        const int& pol,
                        const int& k,
                        const int& m,
                        const int& n);
    static complex_t T3(const arma::cx_mat& field,
                        const int& pol,
                        const int& k,
                        const int& m,
                        const int& n);

    static complex_t T13(const arma::cx_mat& field,
                         const int& pol,
                         const int& k,
                         const int& m,
                         const int& n);

  private:
    int c_radius_;
    int d_radius_;
    bool trained_;

    complex_t intercept_;
    complex_t linear_weight_;
    arma::cx_mat C_;
    arma::cx_mat CC_;
    arma::cx_mat D_;
};
}  // namespace focss

#endif  // FOCSS_EQUALIZER_WDM_SINO_H_
