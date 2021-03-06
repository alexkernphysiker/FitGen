Genetic algorithm for fitting and solving equations (C++17)
===========================================================

This library is distributed under LGPL v.3 license.
Don't be afraid, you can use headers with template classes/functions in your proprietary software if you don't modify them.

What features does this library provide
=======================================

- Solving system of equations (containing any functions)

- Systems of equations containing magnitudes with uncertainties as well

- Fitting points with parametric functions

- You can code you own application of genetic algorithm


Compiling
=========
If you have your git repository with cmake project you can add needed repositories as submodules:

	git submodule add https://github.com/alexkernphysiker/math_h.git
	git submodule add https://github.com/alexkernphysiker/FitGen.git
	git submodule update --init --recursive

Then add to CMakeLists.txt

	add_compile_options(--std=c++17)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread") #option needed for multithreading
	add_definitions(-Dusing_multithread) #option needed for multithreading
	add_subdirectory(math_h)
	add_subdirectory(FitGen)
	include_directories(${MATH_H_INC})
	include_directories(${FITGEN_INC})

Then commit your changes


Examples
========


Solving systems of equations
First example

	Example/system_of_equations.cpp
	
is the example of solving system of equations

	F1(x,y)=G1(x,y)
	F2(x,y)=G2(x,y)

Second example

	Example/system_of_equations2.cpp
	
is the example of solving system of equations containing uncertainties

	F1(x,y)=G1+/-dG1
	F2(x,y)=G2+/-dG2
	
Fitting points with functions
First example

	Example/fit_one_func.cpp
	
is the example of fitting histogram with Gaussian.

Second example

	Example/fit_foreground_background.cpp
	
is the example of fitting points with sum of foreground and background.

