#ifndef LMS_EQUALIZER_H_
#define LMS_EQUALIZER_H_

#include "signal.h"
#include "utility.h"

const double stability_epsilon = 2.2204460492503131e-016;

class LmsEqualizer {
    int order;
    Signal weights;

  public:
    LmsEqualizer();
    LmsEqualizer(const int& filter_order);
    void train(const Signal& desired, const Signal& actual);
    Signal equalize(const Signal& original) const;

  private:
    Signal rotate(const Signal& original, const double& angle) const;
};

#endif  // LMS_EQUALIZER_H_
