cmake_minimum_required(VERSION 3.18)

# set the project name
project(ParTDSE VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

file(GLOB_RECURSE sources src/*.cpp src/include/*.hpp)

# add the executable
add_executable(${PROJECT_NAME} ${sources})

target_compile_options(${PROJECT_NAME} PUBLIC -Wall -pedantic -O3 -march=native -fopenmp)

target_link_options(${PROJECT_NAME} PUBLIC -fopenmp)

target_include_directories(${PROJECT_NAME} PUBLIC src/ src/include/)

set(BLAS_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
    message("OpenBLAS found.")
    include_directories(/opt/OpenBLAS/include/)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${BLAS_LIBRARIES})
endif(BLAS_FOUND)

find_package(yaml-cpp REQUIRED)

find_package(HDF5 REQUIRED)

find_package(GSL REQUIRED)

target_link_libraries(${PROJECT_NAME} PUBLIC yaml-cpp HDF5 GSL)