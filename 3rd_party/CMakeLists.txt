cmake_minimum_required(VERSION 3.20)

project(wigxjpf NONE)

include(FetchContent)

set(
  CMAKE_INSTALL_PREFIX
  "${CMAKE_BINARY_DIR}"
)

set(CMAKE_BUILD_TYPE_INIT "Release")

FetchContent_Declare(${PROJECT_NAME}
    GIT_REPOSITORY    https://github.com/nd-nuclear-theory/wigxjpf.git
    GIT_TAG main
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
)

FetchContent_MakeAvailable(${PROJECT_NAME})