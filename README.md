Genetic algorithm for fitting and solving equations (C++11)
===========================================================
This repository contains the source code of the library, examples and unit tests


Compiling
=========
If you have your git repository with cmake project you can add a submodule

	git submodule add https://github.com/alexkernphysiker/FitGen.git
	git submodule update --init --recursive

Then add to CMakeLists.txt

	add_subdirectory(FitGen)
	include_directories(FitGen/include/)

Then commit your changes :)


CMake Options
=============
	debug
if ON the project is compiled in debug mode

	example
if ON the examples are compiled

	test
if ON the tests are compiled


Examples
========
The examples require gnuplot to be installed for drawing plots of fitting functions.

	Example/equation.cpp
example of using genetic algorithm for solving equation.

	Example/fit_distribution.cpp 
	Example/fit_points.cpp
examples of using genetic algorithm for fitting points with function.
For plotting the results, the examples require gnuplot installed.


Header files
============
	gnuplot.h
	math_h/*.h
links from math_h submodule.

	Genetic/paramset.h
set of parameters.
	
	Genetic/abstract.h
base abstract classes used in the whole library.
	
	Genetic/genetic.h 
template classes that provide different types of genetic algorithm.
	
	Genetic/initialcondition.h 
classes that provide algorithm of generatins points for initialization of population.
	
	Genetic/filter.h 
classes that provide conditions on fitted parameters acting like a filter for new points that appear in the population.
	
	Genetic/equation.h 
template classes and functions providing solving equations using classes declared in other header files.
	
	Genetic/fit.h 
classes needed for fitting algorithms. Contains also classes for plotting 1D fitting results with gnuplot.
	
	Genetic/paramfunc.h 
template classes providing parametric functions for fitting algorithm.


Tests
=====

	tests/*.cpp
The directory tests contains cmake project with unit tests.
The tests require GoogleTest framework to be installed
