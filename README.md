Genetic algorithm for fitting and solving equations (C++11)
===========================================================



Compiling
=========
For using this library, please run:

git clone https://github.com/alexkernphysiker/FitGen.git

cd FitGen

git submodule init

git submodule update

cmake .

make

and you will be able to link the library and run example applications.
For using this library in your git repository please run:

git submodule add https://github.com/alexkernphysiker/FitGen.git

git submodule update --init --recursive

If you have your own cmake project, you can use this library by adding add_subdirectory instruction.
But you will have to add path FitGen/include to your include directories.



CMake Options
=============

debug - if ON the project is compiled in debug mode

example - if ON the examples are compiled

test - if ON the tests are compiled



Examples
========

The examples require gnuplot to be installed for drawing plots of fitting functions.

Example/equation.cpp - example of using genetic algorithm for solving equation.

Example/fit_distribution.cpp and Example/fit_points.cpp - examples of using genetic algorithm for fitting points with function.

For plotting the results, the examples require gnuplot installed


Header files
============

gnuplot.h and math_h/*.h - links from math_h submodule.

Genetic/paramset.h - set of parameters.

Genetic/abstract.h - base abstract classes used in the whole library.

Genetic/genetic.h - template classes that provide different types of genetic algorithm. 
Inheriting one template class from another allows to combine several mutation mechanisms.

Genetic/initialcondition.h - classes that provide algorithm of generatins points for initialization of population.

Genetic/filter.h - classes that provide conditions on fitted parameters acting like a filter for new points that appear in the population.

Genetic/equation.h - template classes and functions providing solving equations using classes declared in other header files.

Genetic/fit.h - classes needed for fitting algorithms. Contains also classes for plotting 1D fitting results with gnuplot.

Genetic/paramfunc.h - template classes providing parametric functions for fitting algorithm. 
Mechanism of inheriting template classes allows to construct complicated ones. 
Parameter number is obtained automatically.

Genetic/paramsort.h - classes for basic features for analysis of parameters distributions. 
Contains classes that allow to split distribution of one parameter into bins and separate ParamSet's by these bins.
More complicated classes allow binning by several parameter indexes.


Tests
=====

The directory tests contains cmake project with unit tests.
The tests require GoogleTest framework to be installed
