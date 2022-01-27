cmake_minimum_required (VERSION 3.16)

project(color_hpp VERSION 0.7.0)

#interface target for CImg
add_library(color_hpp INTERFACE)

#set_target_properties(color_hpp PROPERTIES
#	CXX_STANDARD 17
#	STANDARD_REQUIRED YES)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
	add_subdirectory(example)
endif()
