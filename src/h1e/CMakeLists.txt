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
        ${H5_INC}
)

find_package(yaml-cpp REQUIRED)

find_package(HDF5 REQUIRED)

find_package(OpenMP REQUIRED)

find_package(TBB REQUIRED)

find_package(BLAS REQUIRED HINT ${BLA_HINT})

find_package(GSL REQUIRED)

target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        yaml-cpp hdf5 hdf5_cpp tbb gsl ${BLAS_LIBRARIES} ${HDF5_LIBRARIES}
)
target_link_options(${PROJECT_NAME} PRIVATE -fopenmp)
