cmake_minimum_required(VERSION 3.20)

# set the project name
project(ParTDSE VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

set(SH1_HEADER_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include/dmx_typ.hpp
)

set(SHARED_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cfg_in.cpp
)

set(SHARED_SOURCES2
    ${CMAKE_CURRENT_SOURCE_DIR}/src/fastgl.cpp
)

if(NOT DEFINED BLA_VENDOR)
    set(BLA_VENDOR OpenBLAS)
endif()

if(${BLA_VENDOR} MATCHES OpenBLAS) 
    set(BLAS_INC /opt/OpenBLAS/include)
endif()

if(${BLA_VENDOR} MATCHES Intel10_64lp)
    add_compile_definitions(__INTEL_MKL__)
endif()

if(DEFINED BLA_HINT)
    find_package(BLAS REQUIRED HINT ${BLA_HINT})
    set(BLAS_INC ${BLA_INC})
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")     
    set(Clang_W_FLAG -Wno-c99-extensions) 
endif()

add_subdirectory(3rd_party)

find_package(OpenMP)

add_compile_options(
    -Wall -Wextra -Wshadow -Wnon-virtual-dtor 
    -pedantic -O3 -march=native ${OpenMP_CXX_FLAGS}
    ${Clang_W_FLAG}
)

add_subdirectory(src/h1e)
add_subdirectory(src/w1e)
add_subdirectory(src/d1e)
add_subdirectory(src/tdse)
add_subdirectory(src/pes)
add_subdirectory(src/id2ec)
add_subdirectory(src/d2e)
add_subdirectory(src/r12)
add_subdirectory(src/w2e)
