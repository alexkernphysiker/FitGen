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

Example/fit_distribution.cpp - example of fitting of random numbers distribution.
Here the template classes for parametric function and filtering condition that require std::function parameters are used.

Example/fit_points.cpp - example of fitting data points with x and y errors with complicated function.



Header files
============

paramset.h - set of parameters.

fit_gen.h - base abstract classes used in the whole library.

genetic.h - template classes that provide different types of genetic algorithm. 
Inheriting one template class from another allows to combine several mutation mechanisms.

fitpoints.h - classes representing set of point to fit and several usefull optimality criteria.

paramfunc.h - template classes providing parametric functions. 
Mechanism of inheriting template classes allows to construct complicated ones. 
Parameter number is obtained automatically.
There are also template classes that can accept std::function or other function-like objects.

initialcondition.h - classes that provide algorithm of generatins points for initialization of population.

filter.h - classes that provide conditions on fitted parameters acting like a filter for new points that appear in the population.

equation.h - template classes and functions providing solving equations using classes declared in other header files.
