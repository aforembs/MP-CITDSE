project(h1e)

# add the executable
add_executable(${PROJECT_NAME} ${SHARED_SOURCES2})

target_sources(${PROJECT_NAME}
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/h1e_main.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/h1e.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/bsp_gsl.cpp
)

target_include_directories(${PROJECT_NAME} 
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/../src/include
        ${BLAS_INC}
        ${H5_INC}
)

find_package(HDF5 REQUIRED COMPONENTS C CXX)

find_package(OpenMP REQUIRED)

find_package(TBB REQUIRED)

find_package(BLAS REQUIRED HINT ${BLA_HINT})

find_package(GSL REQUIRED)

target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        yaml-cpp::yaml-cpp tbb ${TBB_LIBRARIES} gsl ${GSL_LIBRARIES} 
        ${BLAS_LIBRARIES} hdf5::hdf5 hdf5::hdf5_cpp
)
target_link_options(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
