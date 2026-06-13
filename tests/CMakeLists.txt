include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.1
)
FetchContent_MakeAvailable(Catch2)

set(TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/test_fastgl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_constants.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_cfg_in.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_bsp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_reference_hydrogen.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_reference_helium.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_reference_tdse.cpp
)

add_executable(tests ${TEST_SOURCES}
    ${CMAKE_SOURCE_DIR}/src/fastgl.cpp
    ${CMAKE_SOURCE_DIR}/src/cfg_in.cpp
)

if(GSL_FOUND)
    target_sources(tests PRIVATE ${CMAKE_SOURCE_DIR}/src/h1e/bsp_gsl.cpp)
endif()

target_include_directories(tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src/include
    ${CMAKE_SOURCE_DIR}/src/h1e/include
    ${CMAKE_BINARY_DIR}/../src/include
    ${H5_INC}
)

find_package(yaml-cpp QUIET)
find_package(HDF5 QUIET COMPONENTS C CXX)
find_package(GSL QUIET)

target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

if(yaml-cpp_FOUND)
    target_link_libraries(tests PRIVATE yaml-cpp::yaml-cpp)
    target_compile_definitions(tests PRIVATE HAS_YAML_CPP)
endif()

if(HDF5_FOUND)
    target_link_libraries(tests PRIVATE hdf5::hdf5 hdf5::hdf5_cpp)
    target_compile_definitions(tests PRIVATE HAS_HDF5)
endif()

if(GSL_FOUND)
    target_link_libraries(tests PRIVATE gsl)
    target_compile_definitions(tests PRIVATE HAS_GSL)
endif()

target_compile_definitions(tests PRIVATE
    TEST_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

target_compile_features(tests PUBLIC cxx_std_17)

include(CTest)
include(Catch)
catch_discover_tests(tests)
