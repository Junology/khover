option(BUILD_DOXYGEN "Generate the documentation using Doxygen." OFF)

if(BUILD_DOXYGEN)
  find_package(Doxygen REQUIRED dot OPTIONAL_COMPONENTS mscgen)

  if(DOXYGEN_FOUND)
	message(${CMAKE_SOURCE_DIR})
	set(DOXY_PROJECT_NAME ${PROJECT_NAME}-doc)
	set(DOXY_README_MD ${CMAKE_SOURCE_DIR}/README.md)
	set(DOXY_PROJECT_INPUT ${CMAKE_SOURCE_DIR}/src)
	set(DOXY_PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/src)

	set(DOXY_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

	configure_file(Doxyfile.in ${DOXY_OUT})

	add_custom_target(doxygen ALL
	  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_OUT}
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	  )
  else(DOXYGEN_FOUND)
	message("Doxygen not found.")
  endif(DOXYGEN_FOUND)
endif(BUILD_DOXYGEN)
