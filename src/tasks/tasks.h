#ifndef TASKS_TASKS_H_
#define TASKS_TASKS_H_

#include <iostream>

void sp_run(const int& n_symbols,
            const double& average_power_dbm,
            std::ostream& set,
            std::ostream& data);

void wdm_run(const int& n_symbols,
             const int& power_index,
             std::string& set,
             std::string& data);

void sp_scoring(const int& n_symbols,
                const double& average_power_dbm,
                std::istream& train,
                std::istream& test,
                std::ostream& data);

void dp_scoring(const int& n_symbols,
                const double& average_power_dbm,
                std::istream& train,
                std::istream& test,
                std::ostream& data);

void wdm_scoring(const int& power_index,
                 const int& span,
                 const std::string& train,
                 const std::string& test,
                 const std::string& data);

void check_model(const int& power_index,
                 const int& span,
                 const std::string& train,
                 const std::string& test,
                 const std::string& data);

void redyuk_scoring();

void find_ber(const int& n_symbols, const int& power_index, std::ostream& data);

void verification(const int& n_samples, const int& n_steps);

#endif  // TASKS_TASKS_H_
