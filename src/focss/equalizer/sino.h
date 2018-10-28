#ifndef FOCSS_EQUALIZER_SINO_H_
#define FOCSS_EQUALIZER_SINO_H_

#include <armadillo>
#include "focss/field.h"

namespace focss {
enum SinoType { ORDINARY, FIXED, HOMOGENEOUS };

class SINO {
  public:
    SINO();
    SINO(const int& radius);
    SINO(const int& radius, const double& lambda);
    SINO(const int& radius, const double& lambda, const SinoType& type);
    SINO& set_radius(const int& radius);
    SINO& set_lambda(const double& lambda);
    SINO& set_type(const SinoType& type);
    SINO& set_linear_weight(const complex_t& linear_weight);

    void train(const Field& transmitted, const Field& received);

    complex_t get_intercept() const;
    complex_t get_linear_weight() const;
    complex_t get_coo_weight() const;
    Field model(const Field& original) const;
    Field equalize(const Field& rx, const Field& tx) const;
    Field linear(const Field& original) const;
    Field cubic(const Field& original) const;

    void save_matrix(const std::string& filename);
    void save_matrix(std::ostream& output);

  private:
    static complex_t term(const Field& field,
                          const int& mode,
                          const int& k,
                          const int& m,
                          const int& n);

  private:
    int radius_;
    double lambda_;
    SinoType type_;
    bool trained_;

    complex_t intercept_;
    complex_t linear_weight_;
    arma::cx_mat weights_;
};
}  // namespace focss

#endif  // FOCSS_EQUALIZER_SINO_H_
