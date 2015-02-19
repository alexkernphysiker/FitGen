cmake_minimum_required(VERSION 2.6)
project(FitGen)
set(CMAKE_CXX_FLAGS "-pthread -std=gnu++11")
add_library(FitGen fit_gen.cpp fitexception.cpp fitpoints.cpp paramrange.cpp params.cpp)
