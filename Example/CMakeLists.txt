#in your own ckame project you should just use 
#add_subdirectory(FitGen)
#for adding library from this repository
cmake_minimum_required(VERSION 2.6)
project(fitapp)
add_definitions(-DUSE_RANDOM_DEVICE)
add_definitions(-DTHREADS_COUNT=2)
include_directories(..)
set(CMAKE_CXX_FLAGS "-pthread -std=gnu++11")
add_executable(fitapp fit.cpp ../fit_gen.cpp ../paramset.cpp ../fitexception.cpp ../fitpoints.cpp ../filter.cpp ../initialconditions.cpp)
set(CMAKE_CXX_FLAGS "-pthread -std=gnu++11")
add_executable(simplefitapp simple_fit.cpp ../fit_gen.cpp ../paramset.cpp ../fitexception.cpp ../fitpoints.cpp ../filter.cpp ../initialconditions.cpp)
set(CMAKE_CXX_FLAGS "-pthread -std=gnu++11")
add_executable(equationapp equation.cpp ../fit_gen.cpp ../paramset.cpp ../fitexception.cpp ../fitpoints.cpp ../filter.cpp ../initialconditions.cpp)
