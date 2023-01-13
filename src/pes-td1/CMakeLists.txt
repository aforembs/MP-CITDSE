project(pes-td1)

# add the executable
add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME}
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/pes1e_main.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/pes1e.cpp
)

target_include_directories(${PROJECT_NAME} 
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

find_package(yaml-cpp REQUIRED)

find_package(HDF5 REQUIRED)

find_package(OpenMP REQUIRED)

target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        yaml-cpp hdf5 hdf5_cpp
)
target_link_options(${PROJECT_NAME} PRIVATE -fopenmp)