project(FitGen)
if(DEFINED MATH_H_INC)
	get_directory_property(has_parent PARENT_DIRECTORY)
	if(has_parent)
		set(FITGEN_INC 
			${CMAKE_CURRENT_SOURCE_DIR}/include
			${MATH_H_INC}
		PARENT_SCOPE)
	endif()

	include_directories(${MATH_H_INC})
	include_directories(include)
	file(GLOB SOURCES src/*.cpp)
	file(GLOB TESTS tests/*.cpp)
	file(GLOB EXAMPLES Example/*.cpp)
	add_library(FitGen SHARED ${SOURCES})

	foreach(example ${EXAMPLES})
		get_filename_component(ex_app ${example} NAME_WE)
		add_executable(${ex_app} ${example})
		target_link_libraries(${ex_app} FitGen)
	endforeach()

	if(DEFINED GTEST)
		add_executable(fitgen_test.exe ${TESTS})
		target_link_libraries(fitgen_test.exe FitGen gtest gtest_main)
	endif()
else()
	message(FATAL_ERROR "FitGen subproject requires math_h subproject to add")
endif()
