project(pes2e)

# add the executable
add_executable(${PROJECT_NAME} ${SHARED_SOURCES} ${SH1_HEADER_FILES})

target_sources(${PROJECT_NAME}
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/pes2e_main.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/pes2e.cpp
)

target_include_directories(${PROJECT_NAME} 
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/../src/include
        ${BLA_INC}
)

find_package(yaml-cpp REQUIRED)

find_package(HDF5 REQUIRED)

find_package(OpenMP REQUIRED)

find_package(TBB REQUIRED)

find_package(BLAS REQUIRED)

target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        yaml-cpp hdf5 hdf5_cpp lapack tbb ${BLAS_LIBRARIES}
)
target_link_options(${PROJECT_NAME} PRIVATE -fopenmp)