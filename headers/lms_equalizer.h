#ifndef LMS_EQUALIZER_H_
#define LMS_EQUALIZER_H_

#include "field.h"
#include "utility.h"

const double stability_epsilon = 2.2204460492503131e-016;

class LmsEqualizer {
    bool trained;
    int length;
    Field weights;

    double mu;
    double alpha;

  public:
    LmsEqualizer();
    LmsEqualizer(const int& filter_length);
    void setFilterLength(const int& filter_length);
    
    void train(const Field& desired, const Field& actual);
    
    Field getWeights() const;
    Field equalize(const Field& original) const;
};

#endif  // LMS_EQUALIZER_H_
