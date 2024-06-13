cmake_minimum_required(VERSION 3.28)
project(bifurcation_analyser)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-std=c++17 -D__USE_FILIB__ -O2 -frounding-math -D__HAVE_MPFR__ -O2 -frounding-math -DFILIB_EXTENDED -DHAVE_SSE -I/usr/local/include)

add_library(libs_module main.cpp
    main.h
    capd_error.h
    neuron.cpp
    neuron.h
    small_funcs.cpp
    small_funcs.h
    solvers.cpp
    solvers.h
)

target_link_libraries(libs_module -L/usr/local/lib)
target_link_libraries(libs_module -lcapd)
target_link_libraries(libs_module -lfilib)
target_link_libraries(libs_module -lmpfr)
target_link_libraries(libs_module -lgmp)

add_executable(bifurcation_analyser main.cpp
    main.h
    capd_error.h
    neuron.cpp
    neuron.h
    small_funcs.cpp
    small_funcs.h
    solvers.cpp
    solvers.h
)

target_link_libraries(bifurcation_analyser libs_module)