#ifndef SINO_H_
#define SINO_H_

#include "focss/field.h"
#include "focss/utility.h"

class SINO {
    unsigned long radius;
    double lambda;
    bool trained;

    Complex intercept;
    ComplexVector weights;

  public:
    SINO();
    SINO(const unsigned long& symbol_radius);
    SINO(const unsigned long& symbol_radius, const double& parameter);
    void setSymbolRadius(const unsigned long& symbol_radius);
    void setRegularizationParameter(const double& parameter);

    void train(const Field& desired, const Field& actual);

    ComplexVector getWeights() const;
    Field equalize(const Field& original) const;

  private:
    Complex term(const Field& field,
                 const unsigned long& center_index,
                 const unsigned long& term_index) const;
};

#endif  // SINO_H_
