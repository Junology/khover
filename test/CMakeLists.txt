cmake_minimum_required(VERSION 3.0)

include_directories(${CMAKE_SOURCE_DIR}/src)

macro(khover_test name srcs)
  add_executable(${name}
	${name}.cpp
	${srcs}
	)
  target_link_libraries(${name} Threads::Threads ${ARGN})

  add_test(
	NAME ${name}
	COMMAND $<TARGET_FILE:${name}>
	)
endmacro(khover_test)

khover_test(
  hnf_test
  ""
  Eigen3::Eigen
  )
