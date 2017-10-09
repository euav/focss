#ifndef LMS_EQUALIZER_H_
#define LMS_EQUALIZER_H_

#include "focss/field.h"
#include "focss/utility.h"

class LeastMeanSquare {
    bool trained;
    unsigned long radius;
    Field weights;

    double mu;
    double alpha;

  public:
    LeastMeanSquare();
    LeastMeanSquare(const unsigned long& symbol_radius);
    void setSymbolRadius(const unsigned long& symbol_radius);
    
    void train(const Field& desired, const Field& actual);
    
    Field getWeights() const;
    Field equalize(const Field& original) const;
};

#endif  // LMS_EQUALIZER_H_
