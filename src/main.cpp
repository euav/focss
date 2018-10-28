#include <fstream>
#include <iostream>
#include <string>
#include "tasks/tasks.h"

int main(int argc, char* argv[]) {
    int sample_size = std::stoi(argv[1]);
    int power_index = std::stoi(argv[2]);

    std::string set = "data/3ch-set-1";
    std::string data = "data/3ch-model-spm";
    std::string train = "data/3ch-set-1";
    std::string test = "data/3ch-set-2";

    std::cout << "Index = " << power_index << std::endl;
    // wdm_run(sample_size, power_index, set, data);
    check_model(power_index, sample_size, train, test, data);

    return 0;
}
