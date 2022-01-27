cmake_minimum_required (VERSION 3.16)

project(color_hpp VERSION 0.7.0)

#interface target for CImg
add_library(color_hpp INTERFACE)
target_include_directories(color_hpp INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_compile_features(color_hpp INTERFACE cxx_std_11)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
	add_subdirectory(example)
endif()
