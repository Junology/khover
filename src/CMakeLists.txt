cmake_minimum_required(VERSION 3.0)

add_executable(
  ${PROJECT_NAME}
  main.cpp
  gausscode.cpp
  linkdiagram.cpp
  cubes.cpp
  khovanov.cpp
)

target_link_libraries(
  ${PROJECT_NAME}
  Threads::Threads
  Eigen3::Eigen
)

set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
