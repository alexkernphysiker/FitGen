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

git submodule init

git submodule update

and add path to it's directory to your include paths.
If you have your own cmake project, you can use this library by adding add_subdirectory instruction.
You can add define USE_RANDOM_DEVICE in the whole your project including this library for using random device instead of pseudorandom number generator.



Examples
========

Example/equation.cpp - example of using genetic algorithm for solving equation.

Example/fit_distribution.cpp and Example/fit_points.cpp - examples of using genetic algorithm for fitting points with function.



Header files
============

paramset.h - set of parameters.

abstract.h - base abstract classes used in the whole library.

genetic.h - template classes that provide different types of genetic algorithm. 
Inheriting one template class from another allows to combine several mutation mechanisms.

initialcondition.h - classes that provide algorithm of generatins points for initialization of population.

filter.h - classes that provide conditions on fitted parameters acting like a filter for new points that appear in the population.

equation.h - template classes and functions providing solving equations using classes declared in other header files.

fit.h - classes needed for fitting algorithms.

paramfunc.h - template classes providing parametric functions for fitting algorithm. 
Mechanism of inheriting template classes allows to construct complicated ones. 
Parameter number is obtained automatically.


Tests
=====

The directory FitGen_tests contains cmake project with unit tests that uses GoogleTest framework.
The GTEST_DIR environment variable must contain path to it's source code (You can download it from official site).
It's subdirectory lib must contain also compiled linkable gtest and gtest_main libraries.
You can compile this project and run the binary to see if you have broken some functions with your changes.
