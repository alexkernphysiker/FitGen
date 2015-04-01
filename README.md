Genetic algorithm for fitting and solving equations (C++11)
===========================================================



Compiling
=========
After clonning this repository, please run:

git submodule init

git submodule update

After that the library can be built with cmake:

cmake .

make

and you will be able to link the library and run example applications.

For using this library in your project, please add the path to cloned repository to your include paths.

If you have your own cmake project, you can use this library by adding add_subdirectory instruction.



Examples
========

Example/equation.cpp - example of using genetic algorithm for solving equation.

Example/fit_distribution.cpp - example of fitting of random numbers distribution.
Here the template classes for parametric function and filtering condition that require std::function parameters are used.

Example/fit_points.cpp - example of fitting data points with x and y errors with complicated function.



Header files
============

paramset.h - set of parameters.

fit_gen.h - base abstract classes used in the whole library.

genetic.h - template classes that provide different types of genetic algorithm. 
Inheriting one template class from another allows to combine several mutation mechanisms.

fitpoints.h - classes that provide optimality criteria for fitting algorithms taking data points into account.
Several template classes allow to extend their functionality.

paramfunc.h - template classes providing parametric functions. 
Mechanism of inheriting template classes allows to construct complicated ones. 
Parameter number is obtained automatically.
There are also template classes that can accept std::function or other function-like objects.

initialcondition.h - classes that provide algorithm of generatins points for initialization of population.

filter.h - classes that provide conditions on fitted parameters acting like a filter for new points that appear in the population.

equation.h - template classes and functions providing solving equations using classes declared in other header files.
