#ifndef LMS_EQUALIZER_H_
#define LMS_EQUALIZER_H_

#include "field.h"
#include "utility.h"

const double stability_epsilon = 2.2204460492503131e-016;

class LeastMeanSquare {
    bool trained;
    unsigned long length;
    Field weights;

    double mu;
    double alpha;

  public:
    LeastMeanSquare();
    LeastMeanSquare(const unsigned long& filter_length);
    void setFilterLength(const unsigned long& filter_length);
    
    void train(const Field& desired, const Field& actual);
    
    Field getWeights() const;
    Field equalize(const Field& original) const;
};

#endif  // LMS_EQUALIZER_H_
