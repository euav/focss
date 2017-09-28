#ifndef SINO_H_
#define SINO_H_

#include "focss/field.h"
#include "focss/utility.h"

class SINO {
    bool trained;
    unsigned long radius;
    Complex scalar_weight;
    Field weights;

  public:
    SINO();
    SINO(const unsigned long& symbol_radius);
    void setSymbolRadius(const unsigned long& symbol_radius);

    void train(const Field& desired, const Field& actual);

    Field getWeights() const;
    Field equalize(const Field& original) const;

  private:
    Complex term(const Field& field,
                 const unsigned long& center_index,
                 const unsigned long& term_index) const;
};

#endif  // SINO_H_
