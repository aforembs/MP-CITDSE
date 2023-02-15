project(r12)

# add the executable
add_executable(${PROJECT_NAME}  
    ${SH1_HEADER_FILES} 
    ${SHARED_SOURCES} 
    ${SHARED_SOURCES2}
)

target_sources(${PROJECT_NAME}
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/r12_main.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/r12.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/integrator.cpp
)

target_include_directories(${PROJECT_NAME} 
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/../src/include
        ${BLAS_INC}
        ${H5_INC}
)

find_package(yaml-cpp REQUIRED)

find_package(HDF5 REQUIRED)

find_package(OpenMP REQUIRED)

find_package(TBB REQUIRED)

find_package(BLAS REQUIRED HINT ${BLA_HINT})

target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        yaml-cpp ${YAML_CPP_LIBRARIES} tbb ${TBB_LIBRARIES} gsl ${GSL_LIBRARIES} 
        wigxjpf ${BLAS_LIBRARIES} hdf5 hdf5_cpp ${HDF5_LIBRARIES}
)
target_link_options(${PROJECT_NAME} PRIVATE -fopenmp)
