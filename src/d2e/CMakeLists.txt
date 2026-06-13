project(d2e)

# add the executable
add_executable(${PROJECT_NAME} ${SHARED_SOURCES})

target_sources(${PROJECT_NAME}
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/dmx2e_main.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/dmx2e.cpp
)

target_include_directories(${PROJECT_NAME} 
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/../src/include
        ${CMAKE_BINARY_DIR}/lib/
        ${H5_INC}
)

find_package(HDF5 REQUIRED COMPONENTS C CXX)

find_package(TBB REQUIRED)

find_package(OpenMP REQUIRED)

target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        yaml-cpp::yaml-cpp tbb 
        ${TBB_LIBRARIES} wigxjpf hdf5::hdf5 hdf5::hdf5_cpp
)
target_link_options(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
