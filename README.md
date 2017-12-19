Genetic algorithm for fitting and solving equations (C++17)
===========================================================

The library is distributed under LGPL v.3 license.
Don't be afraid, you still can use headers containing templates in your proprietary software if you don't modify them.

Question of how LGPL v.3 deals with c++ templates is explained here:

    http://eigen.tuxfamily.org/index.php?title=Licensing_FAQ&oldid=1117

That's a page of C++ library that consists of template headers



Compiling
=========
If you have your git repository with cmake project you can add needed repositories as submodules:

	git submodule add https://github.com/alexkernphysiker/math_h.git
	git submodule add https://github.com/alexkernphysiker/FitGen.git
	git submodule update --init --recursive

Then add to CMakeLists.txt

	add_definitions(--std=c++17)
	add_subdirectory(math_h)
	add_subdirectory(FitGen)
	include_directories(${MATH_H_INC})
	include_directories(${FITGEN_INC})

Then commit your changes

This library still can be compiled with c++11 or c++14 compiler but some features will be less optimized and work slower.



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


Unit tests (require GoogleTest framework)
=========================================

	tests/*.cpp


CMake Options
=============

if ON the tests are compiled

	tests


if ON then the version supporting free threading will be compiled

	threads
