#ifndef TASKS_SYSTEMS_H_
#define TASKS_SYSTEMS_H_

#include "focss/field.h"

void back_to_back(focss::Field& signal);
void cd_compensation(focss::Field& signal);
void cd_compensation(focss::Field& signal, const int& spans);
void forward_propagation(focss::Field& signal);
void forward_propagation(focss::Field& signal, const int& spans);
void backward_propagation(focss::Field& signal, const int& steps = 0);

#endif  // TASKS_SYSTEMS_H_
