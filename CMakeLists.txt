cmake_minimum_required(VERSION 3.3.2)

## Set project
set(PROJECT_NAME "khover")
project(${PROJECT_NAME})
set(PROJECT_VERSION 0.0.1)

set(PROJECT_TARGET_NAME ${PROJECT_NAME})


## Test
enable_testing()


## Language setting
enable_language(CXX)

option(USE_CLANG "Use clang instead of GCC" OFF)

if(USE_CLANG)
  set(CMAKE_C_COMPILER "/usr/bin/clang" CACHE string "clang compiler" FORCE)
  set(CMAKE_CXX_COMPILER "/usr/bin/clang++" CACHE string "clang++ compiler" FORCE)
endif(USE_CLANG)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Debug)


## Libraries
find_package(Threads REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${CMAKE_BINARY_DIR})


## Add subdirectories
add_subdirectory(resources)
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(doc)
