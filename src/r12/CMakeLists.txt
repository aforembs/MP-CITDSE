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

find_package(HDF5 REQUIRED COMPONENTS C CXX)

find_package(OpenMP REQUIRED)

target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        yaml-cpp::yaml-cpp
        wigxjpf hdf5::hdf5 hdf5::hdf5_cpp
)
target_link_options(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
