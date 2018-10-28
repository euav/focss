cmake_minimum_required(VERSION 3.10)
project(FOCSS)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS "-march=native -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -Wall -Wextra -Wshadow -Wno-implicit-fallthrough")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DARMA_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -O3 -DNDEBUG -DARMA_NO_DEBUG")

include_directories(src)
include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_executable(Main
        src/focss/equalizer/admm_regression.cpp
        src/focss/equalizer/phase_shift_equalizer.cpp
        src/focss/equalizer/homothety_equalizer.cpp
        src/focss/equalizer/scalar_equalizer.cpp
        src/focss/equalizer/sino.cpp
        src/focss/equalizer/wdm_sino.cpp
        src/focss/module/amplifier.cpp
        src/focss/module/fiber.cpp
        src/focss/processing/modulation.cpp
        src/focss/solver/ssfm.cpp
        src/focss/definitions.h
        src/focss/field.cpp
        src/focss/functions.cpp
        src/main.cpp
        src/tasks/wdm_scoring.cpp
        src/tasks/dp_scoring.cpp
        src/tasks/sp_scoring.cpp
        src/tasks/sp_run.cpp
        src/tasks/wdm_run.cpp
        src/tasks/redyuk_scoring.cpp
        src/tasks/systems.cpp
        src/tasks/find_ber.cpp src/tasks/verification.cpp src/ppe/ppe.cpp)

target_link_libraries(Main m fftw3 armadillo)
