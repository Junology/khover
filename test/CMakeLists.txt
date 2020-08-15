cmake_minimum_required(VERSION 3.0)

include_directories(${CMAKE_SOURCE_DIR}/src)

macro(khover_test name srcs)
#  message(${srcs})
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
  matrixops_test
  ""
  Eigen3::Eigen
  )

khover_test(
  hnf_test
  ""
  Eigen3::Eigen
  )

khover_test(
  abelian_test
  ""
  Eigen3::Eigen
  )

khover_test(
  hochschild_test
  ""
  Eigen3::Eigen
  )

khover_test(
  states_test
  ""
  Eigen3::Eigen
  )

khover_test(
  gausscode_test
  "${CMAKE_SOURCE_DIR}/src/gausscode.cpp"
  Eigen3::Eigen
  )

khover_test(
  linkdiagram_test
  "${CMAKE_SOURCE_DIR}/src/gausscode.cpp;${CMAKE_SOURCE_DIR}/src/linkdiagram.cpp"
  Eigen3::Eigen
  )
