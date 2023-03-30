project(pes)

# add the executable
add_executable(${PROJECT_NAME} ${SHARED_SOURCES} ${SH1_HEADER_FILES})

target_sources(${PROJECT_NAME}
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/pes_main.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/pes.cpp
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
        yaml-cpp ${YAML_CPP_LIBRARIES} tbb ${TBB_LIBRARIES} 
        ${BLAS_LIBRARIES} hdf5 hdf5_cpp ${HDF5_LIBRARIES}
)
target_link_options(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
