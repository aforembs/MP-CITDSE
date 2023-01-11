project(sort2eEn)

# add the executable
add_executable(${PROJECT_NAME} ${SH1_HEADER_FILES} ${SHARED_SOURCES})

target_sources(${PROJECT_NAME}
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/srt_main.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/srt2e.cpp
)

target_include_directories(${PROJECT_NAME} 
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/../src/include
)

find_package(yaml-cpp REQUIRED)

find_package(HDF5 REQUIRED)

find_package(TBB REQUIRED)

find_package(OpenMP REQUIRED)

target_link_libraries(${PROJECT_NAME} PRIVATE yaml-cpp hdf5 hdf5_cpp tbb)
target_link_options(${PROJECT_NAME} PRIVATE -fopenmp)